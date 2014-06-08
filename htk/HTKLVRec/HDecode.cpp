/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/* developed at:                                               */
/*                                                             */
/*      Machine Intelligence Laboratory                        */
/*      Department of Engineering                              */
/*      University of Cambridge                                */
/*      http://mi.eng.cam.ac.uk/                               */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright:                                          */
/*         2000-2003  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HDecode.c  HTK Large Vocabulary Decoder       */
/* ----------------------------------------------------------- */

char *hdecode_version = "!HVER!HDecode:   3.4.1 [GE 12/03/09]";
char *hdecode_sccs_id = "$Id: HDecode.c,v 1.1.1.1 2006/10/11 09:54:55 jal58 Exp $";

/* this is just the tool that handles command line arguments and
   stuff, all the real magic is in HLVNet and HLVRec */


#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HWave.h"
#include "HLabel.h"
#include "HAudio.h"
#include "HParm.h"
#include "HDict.h"
#include "HModel.h"
#include "HUtil.h"
#include "HTrain.h"
#include "HAdapt.h"
#include "HNet.h"       /* for Lattice */
#include "HLat.h"       /* for Lattice */

#include "config.h"

#include "HLVNet.h"
#include "HLVRec.h"
#include "HLVLM.h"

#include <time.h>

#include <cmdparser.h>

/* -------------------------- Trace Flags & Vars ------------------------ */

#define T_TOP 00001		/* Basic progress reporting */
#define T_OBS 00002		/* Print Observation */
#define T_ADP 00004		/* Adaptation */
#define T_MEM 00010		/* Memory usage, start and finish */

static int trace = 0;

/* -------------------------- Global Variables etc ---------------------- */


static char *langfn;		/* LM filename from commandline */
static char *dictfn;		/* dict filename from commandline */
static char *hmmListfn;		/* model list filename from commandline */
static char *hmmDir = NULL;     /* directory to look for HMM def files */
static char *hmmExt = NULL;     /* HMM def file extension */

static FileFormat ofmt = UNDEFF;	/* Label output file format */
static char *labDir = NULL;	/* output label file directory */
static char *labExt = "rec";	/* output label file extension */
static char *labForm = NULL;	/* output label format */

static Boolean latRescore = FALSE; /* read lattice for each utterance and rescore? */
static char *latInDir = NULL;   /* lattice input directory */
static char *latInExt = "lat";  /* latttice input extension */
static char *latFileMask = NULL; /* mask for reading lattice */

static Boolean latGen = FALSE;  /* output lattice? */
static char *latOutDir = NULL;  /* lattice output directory */
static char *latOutExt = "lat"; /* latttice output extension */
static char *latOutForm = NULL;  /* lattice output format */

static FileFormat dataForm = UNDEFF; /* data input file format */

// static Vocab vocab;		/* wordlist or dictionary */
static HMMSet hset;		/* HMM set */
// static FSLM *lm = NULL;         /* language model */
static LexNet *net = NULL;      /* Lexicon network of all required words/prons */

static char *startWord = "<s>"; /* word used at start of network */
static LabId startLab;          /*   corresponding LabId */
static char *endWord = "</s>";  /* word used at end of network */
static LabId endLab;            /*   corresponding LabId */

static char *spModel = "sp";    /* model used as word end Short Pause */
static LabId spLab;             /*   corresponding LabId */
static char *silModel = "sil";  /* model used as word end Silence */
static LabId silLab;            /*   corresponding LabId */

static Boolean silDict = FALSE; /* does dict contain -/sp/sil variants with probs */

static LogFloat insPen = 0.0;   /* word insertion penalty */

static float acScale = 1.0;     /* acoustic scaling factor */
static float pronScale = 1.0;   /* pronunciation scaling factor */
static float lmScale = 1.0;     /* LM scaling factor */

static int maxModel = 0;        /* max model pruning */
static LogFloat beamWidth = - LZERO;     /* pruning global beam width */
static LogFloat weBeamWidth = - LZERO;   /* pruning wordend beam width */
static LogFloat zsBeamWidth = - LZERO;   /* pruning z-s beam width */
static LogFloat relBeamWidth = - LZERO;  /* pruning relative beam width */
static LogFloat latPruneBeam = - LZERO;  /* lattice pruning beam width */
static LogFloat latPruneAPS = 0;;        /* lattice pruning arcs per sec limit */

static LogFloat fastlmlaBeam = - LZERO;  /* do fast LM la outside this beam */

static int nTok = 32;           /* number of different LMStates per HMM state */
static Boolean useHModel = FALSE; /* use standard HModel OutP functions */
static int outpBlocksize = 1;   /* number of frames for which outP is calculated in one go */
static Observation *obs;        /* array of Observations */

/* transforms/adaptatin */
/* information about transforms */
static XFInfo xfInfo;


/* info for comparing scores from alignment of 1-best with search */
static char *bestAlignMLF;      /* MLF with 1-best alignment */

/* -------------------------- Heaps ------------------------------------- */

// static MemHeap _modelHeap;
// static MemHeap _netHeap;
// static MemHeap _lmHeap;
// static MemHeap _inputBufHeap;
// static MemHeap _transHeap;
// static MemHeap _regHeap;

/* -------------------------- Prototypes -------------------------------- */

void ParseCommandArguments();
void InitAll(int argc, char *argv[]);

void SetConfParms ();
void ReportUsage ();
// DecoderInst *Initialise (void);
// void DoRecognition (DecoderInst *dec, char *fn);
Boolean UpdateSpkrModels (char *fn);

/* linked list storing the info about the 1-best alignment read from BESTALIGNMLF 
   one bestInfo struct per model */
typedef struct _BestInfo BestInfo;
struct _BestInfo {
   int start;           /* frame numbers */
   int end;
   LexNode *ln;
   LLink ll;           /* get rid of this? currently start/end are redundant */
   BestInfo *next;
};

class Vocabulary {
  public:
    Vocabulary() {}

    void init(char* dict_fn) {
      InitVocab (&_vocab);
      if (ReadDict (dict_fn, &_vocab) < SUCCESS)
	HError (9999, "Initialise: ReadDict failed");
    }

    void ConvertSilDict (LabId spLab, LabId silLab, LabId startLab, LabId endLab) {
      ::ConvertSilDict(&_vocab, spLab, silLab, startLab, endLab);
    }

    void MarkAllProns () {
      ::MarkAllProns(&_vocab);
    }

    void MarkAllWords () {
      ::MarkAllWords(&_vocab);
    }

    void UnMarkAllWords () {
      ::UnMarkAllWords(&_vocab);
    }

    void MarkAllWordsfromLat (Lattice* lat, Boolean silDict) {
      ::MarkAllWordsfromLat(&_vocab, lat, silDict);
    }

    Vocab& get_vocab() {
      return _vocab;
    }

  private:
    Vocab _vocab;
};

class LanguageModel {
  public:
    LanguageModel(Vocabulary& vocab): lm(NULL), _vocab(vocab) {
      CreateHeap (&_lmHeap, "LM heap", MSTAK, 1, 0,1000000, 10000000);
    }

    void loadFromFile(char* lm_fn) {
      if (!lm_fn) 
	HError (9999, "HDecode: no LM or lattice specified");
      lm = CreateLM (&_lmHeap, lm_fn, startWord, endWord, &_vocab.get_vocab());
    }

    void loadFromLattice(char* lat_fn, Lattice* lat) {
      lm = CreateLMfromLat (&_lmHeap, lat_fn, lat, &_vocab.get_vocab());
    }

    FSLM* get_lm() {
      return lm;
    }

  private:
    FSLM *lm;         /* language model */
    MemHeap _lmHeap;

    Vocabulary& _vocab;
};

class Decoder {
  public:
    Decoder(LanguageModel& lm, Vocabulary& vocab);
    void init();
    void recognize(char *fn);

    BestInfo *CreateBestInfo (char *fn, HTime frameDur);

    void PrintAlignBestInfo (BestInfo *bestInfo);
    void AnalyseSearchSpace (BestInfo *bestInfo);


    // ===== THIS function does NOT belongs here. Move it somewhere else. =====
    void rescoreLattice(char* fn);

  private:
    LanguageModel& _lm;
    Vocabulary& _vocab;

    DecoderInst* _decInst;

    MemHeap _modelHeap;
    MemHeap _netHeap;
    MemHeap _lmHeap;
    MemHeap _inputBufHeap;
    MemHeap _transHeap;
    MemHeap _regHeap;
};


/* ---------------- Configuration Parameters ---------------------------- */

static ConfParam *cParm[MAXGLOBS];
static int nParm = 0;		/* total num params */


/* ---------------- Debug support  ------------------------------------- */

#if 0
FILE *debug_stdout = stdout;
FILE *debug_stderr = stderr;
#endif

/* ---------------- Process Command Line ------------------------- */

int main (int argc, char *argv[]) {

  Vocabulary vocab;
  LanguageModel languageModel(vocab);
  Decoder decoder(languageModel, vocab);

  InitAll(argc, argv);

  /* load models and initialise decoder */
  vocab.init(dictfn);
  decoder.init();

  /* load 1-best alignment */
  if (bestAlignMLF)
    LoadMasterFile (bestAlignMLF);

  while (NumArgs () > 0) {
    if (NextArg () != STRINGARG)
      HError (4019, "HDecode: Data file name expected");

    char* datafn = GetStrArg ();

    if (trace & T_TOP)
      printf ("File: %s\n", datafn);

    // perform recognition
    decoder.recognize(datafn);
  }

  if (trace & T_MEM) {
    printf ("Memory State on Completion\n");
    PrintAllHeapStats ();
  }

  /* maybe output transforms for last speaker */
  UpdateSpkrStats(&hset, &xfInfo, NULL); 

  Exit(0);             /* maybe print config and exit */
  return 0;
}

/* SetConfParms: set conf parms relevant to this tool */
void SetConfParms () {
   int i;
   double f;
   Boolean b;
   char buf[MAXSTRLEN];

   nParm = GetConfig ("HDECODE", TRUE, cParm, MAXGLOBS);
   if (nParm > 0) {
      if (GetConfInt (cParm, nParm, "TRACE", &i))
	 trace = i;
      if (GetConfStr (cParm, nParm, "STARTWORD", buf))
         startWord = CopyString (&gstack, buf);
      if (GetConfStr (cParm, nParm, "ENDWORD", buf))
         endWord = CopyString (&gstack, buf);
      if (GetConfFlt (cParm, nParm, "LATPRUNEBEAM", &f))
         latPruneBeam  = f;
      if (GetConfFlt (cParm, nParm, "FASTLMLABEAM", &f))
         fastlmlaBeam  = f;
      if (GetConfFlt (cParm, nParm, "LATPRUNEAPS", &f))
         latPruneAPS  = f;
      if (GetConfStr (cParm, nParm, "BESTALIGNMLF", buf))
         bestAlignMLF = CopyString (&gstack, buf);
      if (GetConfBool (cParm, nParm, "USEHMODEL",&b)) useHModel = b;
      if (GetConfStr(cParm,nParm,"LATFILEMASK",buf)) {
         latFileMask = CopyString(&gstack, buf);
      }
   }
}

void ReportUsage () {
   printf ("\nUSAGE: HDecode [options] VocabFile HMMList DataFiles...\n\n");
   printf (" Option                                   Default\n\n");
   printf (" -m      enable XForm and use inXForm        off\n");

   printf (" -d s    dir to find hmm definitions       current\n");
   printf (" -i s    Output transcriptions to MLF s      off\n");
   printf (" -k i    block size for outP calculation     1\n");
   printf (" -l s    dir to store label files	    current\n");
   printf (" -o s    output label formating NCSTWMX      none\n");
   printf (" -h s    speaker name pattern                none\n");
   printf (" -p f    word insertion penalty              0.0\n");
   printf (" -a f    acoustic scale factor               1.0\n");
   printf (" -r f    pronunciation scale factor          1.0\n");
   printf (" -s f    LM scale factor                     1.0\n");
   printf (" -t f    pruning beam width                  none\n");
   printf (" -u i    max model pruning                   0\n");
   printf (" -v f    wordend beam width                  0.0\n");
   printf (" -n i    number of tokens per state          32\n");
   printf (" -w s    use language model                  none\n");
   printf (" -x s    extension for hmm files             none\n");
   printf (" -y s    output label file extension         rec\n");
   printf (" -z s    generate lattices with extension s  off\n");
   printf (" -q s    output lattices format ABtvaldmnr  tvaldmr\n");
   printf (" -R s    best align MLF                      off\n");
   printf (" -X ext  set input lattice extension         lat\n");
   PrintStdOpts ("EJFHLSTP");
   printf ("\n\n");

   printf ("build-time options: ");
#ifdef MODALIGN
   printf ("MODALIGN ");
#endif   
#ifdef TSIDOPT
   printf ("TSIDOPT ");
#endif   
   printf ("\n  sizes: PronId=%d  LMId=%d \n", sizeof (PronId), sizeof (LMId));
}


// DecoderInst *Initialise (void) { }


/**********  align best code  ****************************************/

/* find the LN_MODEL lexnode following ln that has label lab
   step over LN_CON and LN_WORDEND nodes.
   return NULL if not found
*/
BestInfo *FindLexNetLab (MemHeap *heap, LexNode *ln, LLink ll, HTime frameDur)
{
   int i;
   LexNode *follLN;
   MLink m;
   BestInfo *info, *next;

   if (!ll->succ) {
      info = (BestInfo*) New (heap, sizeof (BestInfo));
      info->next = NULL;
      info->ll = NULL;
      info->ln = NULL;
      info->start = info->end = 0;
      return info;
   }
   
   for (i = 0; i < ln->nfoll; ++i) {
      follLN = ln->foll[i];
      if (follLN->type == LN_MODEL) {
         m = FindMacroStruct (&hset, 'h', follLN->data.hmm);
         if (m->id == ll->labid) {
            /*            printf ("found  %8.0f %8.0f %8s  %p\n", ll->start, ll->end, ll->labid->name, follLN); */
            next = FindLexNetLab (heap, follLN, ll->succ, frameDur);
            if (next) {
               info = (BestInfo*) New (heap, sizeof (BestInfo));
               info->next = next;
               info->start = ll->start / (frameDur*1.0e7);
               info->end = ll->end / (frameDur*1.0e7);
               info->ll = ll;
               info->ln = follLN;
               return info;
            }
            /*            printf ("damn got 0 back searching for %8s\n", ll->labid->name); */
         }
      }
      else {
         /*         printf ("searching for %8s recursing\n", ll->labid->name); */
         next = FindLexNetLab (heap, follLN, ll, frameDur);
         if (next) {
            info = (BestInfo*) New (heap, sizeof (BestInfo));
            info->next = next;
            info->start = info->end = ll->start / (frameDur*1.0e7);
            info->ll = ll;
            info->ln = follLN;
            return info;
         }
         /*         printf ("damn got 0 back from recursion\n"); */
      }
   }
   
   return NULL;
}

void InitAll(int argc, char *argv[]) {

  if (InitShell (argc, argv, hdecode_version, hdecode_sccs_id) < SUCCESS)
    HError (4000, "HDecode: InitShell failed");

  InitMem ();
  InitMath ();
  InitSigP ();
  InitWave ();
  InitLabel ();
  InitAudio ();
  InitModel ();
  if (InitParm () < SUCCESS)
    HError (4000, "HDecode: InitParm failed");
  InitUtil ();
  InitDict ();
  InitLVNet ();
  InitLVLM ();
  InitLVRec ();
  InitAdapt (&xfInfo);
  InitLat ();

  if (!InfoPrinted () && NumArgs () == 0) ReportUsage ();
  if (NumArgs () == 0) Exit (0);

  SetConfParms ();
  ParseCommandArguments();
}

void ParseCommandArguments() {
   while (NextArg () == SWITCHARG) {
      char *s = GetSwtArg ();
      if (strlen (s) != 1)
	 HError (4019, "HDecode: Bad switch %s; must be single letter", s);
      switch (s[0]) {
      case 'd':
	 if (NextArg() != STRINGARG)
	    HError(4119,"HDecode: HMM definition directory expected");
	 hmmDir = GetStrArg(); 
	 break;
      case 'x':
	 if (NextArg() != STRINGARG)
	    HError(4119,"HDecode: HMM file extension expected");
	 hmmExt = GetStrArg(); 
	 break;
	 
      case 'i':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output MLF file name expected");
	 if (SaveToMasterfile (GetStrArg ()) < SUCCESS)
	    HError (4014, "HDecode: Cannot write to MLF");
	 break;

      case 'P':
	 if (NextArg () != STRINGARG)
	    HError (3219, "HVite: Target Label File format expected");
	 if ((ofmt = Str2Format (GetStrArg ())) == ALIEN)
	    HError (-3289,
		    "HVite: Warning ALIEN Label output file format set");
	 break;

      case 'l':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Label/Lattice output directory expected");
	 labDir = GetStrArg ();
         latOutDir = labDir;
	 break;
      case 'o':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output label format expected");
	 labForm = GetStrArg ();
	 break;
      case 'y':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output label file extension expected");
	 labExt = GetStrArg ();
	 break;

      case 'X':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Input Lattice file extension expected");
	 latInExt = GetStrArg ();
	 break;
      case 'L':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Input Lattice directory expected");
	 latInDir = GetStrArg ();
	 break;

      case 'q':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output lattice format expected");
	 latOutForm = GetStrArg ();
	 break;
      case 'z':
         latGen = TRUE;
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output lattice file extension expected");
	 latOutExt = GetStrArg ();
	 break;

      case 'p':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: word insertion penalty expected");
         insPen = GetFltArg ();
         if (insPen > 0.0)
            HError (-1, "HDecode: positive word insertion penalty???");
	 break;
      case 'a':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: acoustic scale factor expected");
	 acScale = GetFltArg ();
	 break;
      case 'r':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: pronunciation scale factor expected");
	  pronScale = GetFltArg ();
          silDict = TRUE;       /* #### maybe separate switch for this? */
	 break;
      case 's':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: LM scale factor expected");
         lmScale= GetFltArg ();
	 break;


      case 'u':
	 if (NextArg () != INTARG)
	    HError (4019, "HDecode: max model pruning limit expected");
         maxModel = GetIntArg ();
	 break;

      case 't':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: beam width expected");
	 beamWidth = GetFltArg ();
         if (latPruneBeam == -LZERO)
            latPruneBeam = beamWidth;
         relBeamWidth = beamWidth;
         if (NextArg () == FLOATARG)
            relBeamWidth = GetFltArg ();
	 break;

      case 'v':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: wordend beam width expected");
         weBeamWidth = GetFltArg ();
         zsBeamWidth = weBeamWidth;
	 if (NextArg () == FLOATARG)
            zsBeamWidth = GetFltArg ();
         break;

      case 'w':
	 if (NextArg() != STRINGARG) {
            /*	    HError (4119,"HDecode: LM File name expected"); */
            latRescore = TRUE;
         }
         else
            langfn = GetStrArg();
	 break;

      case 'n':
	 nTok = GetChkedInt (0, 1024, s);
	 break;

      case 'k':
	 outpBlocksize = GetChkedInt (0, MAXBLOCKOBS, s);
	 break;

      case 'H':
	 if (NextArg() != STRINGARG)
	    HError (4119,"HDecode: MMF File name expected");
	 AddMMF (&hset, GetStrArg()); 
	 break;
      case 'T':
	 trace = GetChkedInt (0, 1000, s);
	 break;

      case 'h':
         if (NextArg()!=STRINGARG)
	    HError (4019, "HDecode: Speaker name pattern expected");
         xfInfo.outSpkrPat = GetStrArg();
         if (NextArg()==STRINGARG) {
            xfInfo.inSpkrPat = GetStrArg();
            if (NextArg()==STRINGARG)
               xfInfo.paSpkrPat = GetStrArg(); 
         }
         if (NextArg() != SWITCHARG)
	    HError (4019, "HDecode: cannot have -h as the last option");
         break;
      case 'm':
	 xfInfo.useInXForm = TRUE;
         break;
      case 'E':
         if (NextArg()!=STRINGARG)
            HError(4019,"HDecode: parent transform directory expected");
	 xfInfo.usePaXForm = TRUE;
         xfInfo.paXFormDir = GetStrArg(); 
         if (NextArg()==STRINGARG)
            xfInfo.paXFormExt = GetStrArg(); 
	 if (NextArg() != SWITCHARG)
            HError(4019,"HDecode: cannot have -E as the last option");	  
         break;              
      case 'J':
         if (NextArg()!=STRINGARG)
            HError(4019,"HDecode: input transform directory expected");
         AddInXFormDir(&hset,GetStrArg());
         if (NextArg()==STRINGARG)
            xfInfo.inXFormExt = GetStrArg(); 
	 if (NextArg() != SWITCHARG)
            HError(4019,"HDecode: cannot have -J as the last option");	  
         break;              
      case 'K':
         HError(4019,"HDecode: transform estimation (-K option) not supported yet");	  
         if (NextArg()!=STRINGARG)
            HError(4019,"HDecode: output transform directory expected");
         xfInfo.outXFormDir = GetStrArg(); 
	 xfInfo.useOutXForm = TRUE;
         if (NextArg()==STRINGARG)
            xfInfo.outXFormExt = GetStrArg(); 
	 if (NextArg() != SWITCHARG)
            HError(4019,"HDecode: cannot have -K as the last option");	  
         break;              
      case 'N':
         HError (4019, "HDecode: old style fv transform not supported!");
	 break;
      case 'Q':
         HError (4019, "HDecode: old style mllr transform not supported!");
	 break;

      case 'R':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: best align MLF name expected");
	 bestAlignMLF = GetStrArg ();
	 break;


      default:
	 HError (4019, "HDecode: Unknown switch %s", s);
      }
   }

   if (NextArg () != STRINGARG)
      HError (4019, "HDecode Vocab file name expected");
   dictfn = GetStrArg ();

   if (NextArg () != STRINGARG)
      HError (4019, "HDecode model list file name expected");
   hmmListfn = GetStrArg ();

   if (beamWidth > -LSMALL)
      HError (4019, "main beam is too wide!");

   if (xfInfo.useInXForm) {
      if (!useHModel) {
         HError (-4019, "HDecode: setting USEHMODEL to TRUE.");
         useHModel = TRUE;
      }
      if (outpBlocksize != 1) {
         HError (-4019, "HDecode: outP blocksize >1 not supported with new XForm code! setting to 1.");
         outpBlocksize = 1;
      }
   }   
}

/*****************  main recognition function  ************************/

// void DoRecognition (DecoderInst *dec, char *fn) { }

Boolean UpdateSpkrModels (char *fn) {
   HError (1, "MLLR or FV transforms not supported");
   return FALSE;
}

// =======================================================

Decoder::Decoder(LanguageModel& lm, Vocabulary& vocab): _lm(lm), _vocab(vocab), _decInst(NULL) {
   /* init model heap & set early to support loading MMFs */
   CreateHeap(&_modelHeap, "Model heap",  MSTAK, 1, 0.0, 100000, 800000 );
   CreateHMMSet(&hset, &_modelHeap, TRUE); 
}

void Decoder::init() {

   Boolean eSep;
   Boolean modAlign;

   /* init Heaps */
   CreateHeap (&_netHeap, "Net heap", MSTAK, 1, 0,100000, 800000);
   CreateHeap (&_lmHeap, "LM heap", MSTAK, 1, 0,1000000, 10000000);
   CreateHeap (&_transHeap,"Transcription heap",MSTAK,1,0,8000,80000);

   /* Read dictionary */
   if (trace & T_TOP)
      printf ("Reading dictionary from %s\n", dictfn);

   /* Read accoustic models */
   if (trace & T_TOP)
      printf ("Reading acoustic models...");

   if (MakeHMMSet (&hset, hmmListfn) < SUCCESS) 
      HError (4128, "Initialise: MakeHMMSet failed");
   if (LoadHMMSet (&hset, hmmDir, hmmExt) < SUCCESS) 
      HError (4128, "Initialise: LoadHMMSet failed");
   
   /* convert to INVDIAGC */
   ConvDiagC (&hset, TRUE);
   ConvLogWt (&hset);
   
   if (trace&T_TOP)
      printf("Read %d physical / %d logical HMMs\n", hset.numPhyHMM, hset.numLogHMM);  

   /* process dictionary */
   startLab = GetLabId (startWord, FALSE);
   if (!startLab) 
      HError (9999, "HDecode: cannot find STARTWORD '%s'\n", startWord);
   endLab = GetLabId (endWord, FALSE);
   if (!endLab) 
      HError (9999, "HDecode: cannot find ENDWORD '%s'\n", endWord);

   spLab = GetLabId (spModel, FALSE);
   if (!spLab)
      HError (9999, "HDecode: cannot find label 'sp'");
   silLab = GetLabId (silModel, FALSE);
   if (!silLab)
      HError (9999, "HDecode: cannot find label 'sil'");

   if (silDict) {    /* dict contains -/sp/sil variants (with probs) */
      _vocab.ConvertSilDict (spLab, silLab, startLab, endLab);

      /* check for skip in sp model */
      { 
         LabId spLab;
         HLink spHMM;
         MLink spML;
         int N;

         spLab = GetLabId ("sp", FALSE);
         if (!spLab)
            HError (9999, "cannot find 'sp' model.");

         spML = FindMacroName (&hset, 'l', spLab);
         if (!spML)
            HError (9999, "cannot find model for sp");
         spHMM = (HLink) spML->structure;
         N = spHMM->numStates;

         if (spHMM->transP[1][N] > LSMALL)
            HError (9999, "HDecode: using -/sp/sil dictionary but sp contains tee transition!");
      }
   }
   else {       /* lvx-style dict (no sp/sil at wordend */
     _vocab.MarkAllProns();
   }

   if (!latRescore) {

      /* mark all words  for inclusion in Net */
     _vocab.MarkAllWords();

      /* create network */
      net = CreateLexNet (&_netHeap, &_vocab.get_vocab(), &hset, startWord, endWord, silDict);
      
      /* Read language model */
      if (trace & T_TOP)
         printf ("Reading language model from %s\n", langfn);
      
      _lm.loadFromFile(langfn);
   }

   modAlign = FALSE;
   if (latOutForm) {
      if (strchr (latOutForm, 'd'))
         modAlign = TRUE;
      if (strchr (latOutForm, 'n'))
         HError (9999, "DoRecognition: likelihoods for model alignment not supported");
   }

   /* create Decoder instance */
   _decInst = CreateDecoderInst (&hset, _lm.get_lm(), nTok, TRUE, useHModel, outpBlocksize,
                            bestAlignMLF ? TRUE : FALSE,
                            modAlign);
   
   /* create buffers for observations */
   SetStreamWidths (hset.pkind, hset.vecSize, hset.swidth, &eSep);

   obs = (Observation *) New (&gcheap, outpBlocksize * sizeof (Observation));
   for (int i = 0; i < outpBlocksize; ++i)
      obs[i] = MakeObservation (&gcheap, hset.swidth, hset.pkind, 
                                (hset.hsKind == DISCRETEHS), eSep);

   CreateHeap (&_inputBufHeap, "Input Buffer Heap", MSTAK, 1, 1.0, 80000, 800000);

   /* Initialise adaptation */

   /* sort out masks just in case using adaptation */
   if (xfInfo.inSpkrPat == NULL) xfInfo.inSpkrPat = xfInfo.outSpkrPat; 
   if (xfInfo.paSpkrPat == NULL) xfInfo.paSpkrPat = xfInfo.outSpkrPat; 

   if (xfInfo.useOutXForm) {
      CreateHeap(&_regHeap,   "regClassStore",  MSTAK, 1, 0.5, 1000, 8000 );
      /* This initialises things - temporary hack - THINK!! */
      CreateAdaptXForm(&hset, "tmp");

      /* online adaptation not supported yet! */
   }
}

void Decoder::rescoreLattice(char* fn) {
  /* read lattice and create LM */
  char latfn[MAXSTRLEN],buf3[MAXSTRLEN];
  FILE *latF;
  Boolean isPipe;
  Lattice *lat;

  /* clear out previous LexNet, Lattice and LM structures */
  ResetHeap (&_lmHeap);
  ResetHeap (&_netHeap);

  if (latFileMask != NULL ) { /* support for rescoring lattoce masks */
    if (!MaskMatch(latFileMask, buf3 , fn))
      HError(2319,"HDecode: mask %s has no match with segemnt %s",latFileMask,fn);
    MakeFN (buf3, latInDir, latInExt, latfn);
  } else {
    MakeFN (fn, latInDir, latInExt, latfn);
  }

  if (trace & T_TOP)
    printf ("Loading Lattice from %s\n", latfn);
  {
    latF = FOpen (latfn, NetFilter, &isPipe);
    if (!latF)
      HError (9999, "DoRecognition: Cannot open lattice file %s\n", latfn);

    /* #### maybe separate lattice heap? */
    lat = ReadLattice (latF, &_lmHeap, &_vocab.get_vocab(), FALSE, FALSE);
    FClose (latF, isPipe);
    if (!lat)
      HError (9999, "DoRecognition: cannot read lattice file %s\n", latfn);
  }

  /* mark prons of all words in lattice */
  _vocab.UnMarkAllWords();
  _vocab.MarkAllWordsfromLat(lat, silDict);

  /* create network of all the words/prons marked (word->aux and pron->aux == 1) */
  if (trace & T_TOP)
    printf ("Creating network\n");
  net = CreateLexNet (&_netHeap, &_vocab.get_vocab(), &hset, startWord, endWord, silDict);

  /* create LM based on pronIds defined by CreateLexNet */
  if (trace & T_TOP)
    printf ("Creating language model\n");

  _lm.loadFromLattice(latfn, lat);
  _decInst->lm = _lm.get_lm();
}

void Decoder::recognize(char *fn) {

   char buf1[MAXSTRLEN], buf2[MAXSTRLEN];
   ParmBuf parmBuf;
   BufferInfo pbInfo;
   int frameN, frameProc, i, bs;
   Transcription *trans;
   Lattice *lat;
   clock_t startClock, endClock;
   double cpuSec;
   Observation *obsBlock[MAXBLOCKOBS];
   BestInfo *bestAlignInfo = NULL;

   /* This handles the initial input transform, parent transform setting
      and output transform creation */
   { 
      Boolean changed;
      changed = UpdateSpkrStats(&hset, &xfInfo, fn);
   }

   startClock = clock();

   /* get transcrition of 1-best alignment */
   if (bestAlignMLF)
      bestAlignInfo = CreateBestInfo (fn, pbInfo.tgtSampRate/1.0e7);
   
   parmBuf = OpenBuffer (&_inputBufHeap, fn, 50, dataForm, TRI_UNDEF, TRI_UNDEF);
   if (!parmBuf)
      HError (9999, "HDecode: Opening input failed");
   
   GetBufferInfo (parmBuf, &pbInfo);
   if (pbInfo.tgtPK != hset.pkind)
      HError (9999, "HDecode: Incompatible parm kinds %s vs. %s",
              ParmKind2Str (pbInfo.tgtPK, buf1),
              ParmKind2Str (hset.pkind, buf2));
              
   if (latRescore)
     rescoreLattice(fn);

   if (weBeamWidth > beamWidth)
      weBeamWidth = beamWidth;
   if (zsBeamWidth > beamWidth)
      zsBeamWidth = beamWidth;

   InitDecoderInst (_decInst, net, pbInfo.tgtSampRate, beamWidth, relBeamWidth,
                    weBeamWidth, zsBeamWidth, maxModel,
                    insPen, acScale, pronScale, lmScale, fastlmlaBeam);

   net->vocabFN = dictfn;
   _decInst->utterFN = fn;

   frameN = frameProc = 0;
   while (BufferStatus (parmBuf) != PB_CLEARED) {
      ReadAsBuffer (parmBuf, &obs[frameN % outpBlocksize]);

      if (frameN+1 >= outpBlocksize) {  /* enough frames available */
         if (trace & T_OBS)
            PrintObservation (frameProc+1, &obs[frameProc % outpBlocksize], 13);
         for (i = 0; i < outpBlocksize; ++i)
            obsBlock[i] = &obs[(frameProc + i) % outpBlocksize];

#ifdef DEBUG_TRACE
         fprintf(stdout, "\nProcessing frame %d :\n", frameProc);
#endif
         
         ProcessFrame (_decInst, obsBlock, outpBlocksize, xfInfo.inXForm);
         if (bestAlignInfo)
            this->AnalyseSearchSpace (bestAlignInfo);
         ++frameProc;
      }
      ++frameN;
   }
   CloseBuffer (parmBuf);

   /* process remaining frames (no full blocks available anymore) */
   for (bs = outpBlocksize-1; bs >=1; --bs) {
      if (trace & T_OBS)
         PrintObservation (frameProc+1, &obs[frameProc % outpBlocksize], 13);
      for (i = 0; i < bs; ++i)
         obsBlock[i] = &obs[(frameProc + i) % outpBlocksize];
      
      ProcessFrame (_decInst, obsBlock, bs, xfInfo.inXForm);
      if (bestAlignInfo)
         this->AnalyseSearchSpace (bestAlignInfo);
      ++frameProc;
   }
   assert (frameProc == frameN);

   
   endClock = clock();
   cpuSec = (endClock - startClock) / (double) CLOCKS_PER_SEC;
   printf ("CPU time %f  utterance length %f  RT factor %f\n",
           cpuSec, frameN*_decInst->frameDur, cpuSec / (frameN*_decInst->frameDur));

   trans = TraceBack (&_transHeap, _decInst);

   /* save 1-best transcription */
   /* the following is from HVite.c */
   if (trans) {
      char labfn[MAXSTRLEN];

      if (labForm != NULL)
         ReFormatTranscription (trans, pbInfo.tgtSampRate, FALSE, FALSE,
                                strchr(labForm,'X')!=NULL,
                                strchr(labForm,'N')!=NULL,strchr(labForm,'S')!=NULL,
                                strchr(labForm,'C')!=NULL,strchr(labForm,'T')!=NULL,
                                strchr(labForm,'W')!=NULL,strchr(labForm,'M')!=NULL);
      
      MakeFN (fn, labDir, labExt, labfn);

      if (LSave (labfn, trans, ofmt) < SUCCESS)
         HError(9999, "DoRecognition: Cannot save file %s", labfn);
      if (trace & T_TOP)
         PrintTranscription (trans, "1-best hypothesis");

      Dispose (&_transHeap, trans);
   }

   if (latGen) {
      lat = LatTraceBack (&_transHeap, _decInst);

      /* prune lattice */
      if (lat && latPruneBeam < - LSMALL) {
         lat = LatPrune (&_transHeap, lat, latPruneBeam, latPruneAPS);
      }

      /* the following is from HVite.c */
      if (lat) {
         char latfn[MAXSTRLEN];
         char *p;
         Boolean isPipe;
         FILE *file;
         LatFormat form;
         
         MakeFN (fn, latOutDir, latOutExt, latfn);
         file = FOpen (latfn, NetOFilter, &isPipe);
         if (!file) 
            HError (999, "DoRecognition: Could not open file %s for lattice output",latfn);
         if (!latOutForm)
            form = (HLAT_DEFAULT & ~HLAT_ALLIKE)|HLAT_PRLIKE;
         else {
            for (p = latOutForm, form=0; *p != 0; p++) {
               switch (*p) {
               case 'A': form|=HLAT_ALABS; break;
               case 'B': form|=HLAT_LBIN; break;
               case 't': form|=HLAT_TIMES; break;
               case 'v': form|=HLAT_PRON; break;
               case 'a': form|=HLAT_ACLIKE; break;
               case 'l': form|=HLAT_LMLIKE; break;
               case 'd': form|=HLAT_ALIGN; break;
               case 'm': form|=HLAT_ALDUR; break;
               case 'n': form|=HLAT_ALLIKE; 
                  HError (9999, "DoRecognition: likelihoods for model alignment not supported");
                  break;
               case 'r': form|=HLAT_PRLIKE; break;
               }
            }
         }
         if (WriteLattice (lat, file, form) < SUCCESS)
            HError(9999, "DoRecognition: WriteLattice failed");
         
         FClose (file,isPipe);
         Dispose (&_transHeap, lat);
      }
   }


#ifdef COLLECT_STATS
   printf ("Stats: nTokSet %lu\n", _decInst->stats.nTokSet);
   printf ("Stats: TokPerSet %f\n", _decInst->stats.sumTokPerTS / (double) _decInst->stats.nTokSet);
   printf ("Stats: activePerFrame %f\n", _decInst->stats.nActive / (double) _decInst->stats.nFrames);
   printf ("Stats: activateNodePerFrame %f\n", _decInst->stats.nActivate / (double) _decInst->stats.nFrames);
   printf ("Stats: deActivateNodePerFrame %f\n\n", 
           _decInst->stats.nDeActivate / (double) _decInst->stats.nFrames);
#if 0
   printf ("Stats: LMlaCacheHits %ld\n", _decInst->stats.nLMlaCacheHit);
   printf ("Stats: LMlaCacheMiss %ld\n", _decInst->stats.nLMlaCacheMiss);
#endif
#ifdef COLLECT_STATS_ACTIVATION
   {
      int i;
      for (i = 0; i <= STATS_MAXT; ++i)
         printf ("T %d Dead %lu Live %lu\n", i, _decInst->stats.lnDeadT[i], _decInst->stats.lnLiveT[i]);
   }
#endif
#endif

   if (trace & T_MEM) {
      printf ("memory stats at end of recognition\n");
      PrintAllHeapStats ();
   }

   ResetHeap (&_inputBufHeap);
   ResetHeap (&_transHeap);
   CleanDecoderInst (_decInst);
}

BestInfo* Decoder::CreateBestInfo (char *fn, HTime frameDur)
{
   char alignFN[MAXFNAMELEN];
   Transcription *bestTrans;
   LLink ll;
   LexNode *ln;
   MLink m;
   LabId lnLabId;
   BestInfo *bestAlignInfo;

   MakeFN (fn, "", "rec", alignFN);
   bestTrans = LOpen (&_transHeap, alignFN, HTK);
      
   /* delete 'sp' or 'sil' before final 'sil' if it is there
      these are always inserted by HVite but not possible in HDecode's net structure*/
   if (bestTrans->head->tail->pred->pred->labid == spLab ||
       bestTrans->head->tail->pred->pred->labid == silLab) {
      LLink delLL;
      
      delLL = bestTrans->head->tail->pred->pred;
      /* add sp's frames (if any) to final sil */
      delLL->succ->start = delLL->pred->end;
      
      delLL->pred->succ = delLL->succ;
      delLL->succ->pred = delLL->pred;
   }
   
   ln = net->start;
   assert (ln->type == LN_MODEL);
   m = FindMacroStruct (&hset, 'h', ln->data.hmm);
   lnLabId = m->id;
   
   /* info for net start node */
   ll = bestTrans->head->head->succ;

   assert (ll->labid == lnLabId);
   bestAlignInfo = (BestInfo*) New (&_transHeap, sizeof (BestInfo));
   bestAlignInfo->start = ll->start / (frameDur*1.0e7);
   bestAlignInfo->end = ll->end / (frameDur*1.0e7);
   bestAlignInfo->ll = ll;
   bestAlignInfo->ln = ln;
   
   /* info for all the following nodes */
   bestAlignInfo->next = FindLexNetLab (&_transHeap, ln, ll->succ, frameDur);
   
   {
      BestInfo *b;
      for (b = bestAlignInfo; b->next; b = b->next)
         printf ("%d %d %8s %p\n", b->start, b->end, b->ll->labid->name, b->ln);
   }

   return bestAlignInfo;
}

void Decoder::PrintAlignBestInfo (BestInfo *b)
{
   LexNodeInst *inst;
   TokScore score;
   LabId monoPhone;
   LogDouble phonePost;

   inst = b->ln->inst;
   score = inst ? inst->best : LZERO;

   if (b->ln->type == LN_MODEL) {
      monoPhone =(LabId) b->ln->data.hmm->hook;
      phonePost = _decInst->phonePost[(int) monoPhone->aux];
   } else
      phonePost = 999.99;

   int l = _decInst->nLayers-1;
   while (_decInst->net->layerStart[l] > b->ln) {
      --l;
      assert (l >= 0);
   }
   
   printf ("BESTALIGN frame %4d best %.3f alignbest %d -> %d ln %p layer %d score %.3f phonePost %.3f\n", 
           _decInst->frame, _decInst->bestScore, 
           b->start, b->end, b->ln, l, score, phonePost);
}

void Decoder::AnalyseSearchSpace (BestInfo *bestInfo) {
  if (!bestInfo)
    return;

  BestInfo *b;
  LabId monoPhone;

  monoPhone =(LabId) _decInst->bestInst->node->data.hmm->hook;
  printf ("frame %4d best %.3f phonePost %.3f\n", _decInst->frame, 
      _decInst->bestScore, _decInst->phonePost[(int) monoPhone->aux]);

  for (b = bestInfo; b; b = b->next) {
    if (b->start < _decInst->frame && b->end >= _decInst->frame) 
      break;
  }
  if (b) {
    this->PrintAlignBestInfo (b);
    for (b = b->next; b && b->start == b->end && b->start == _decInst->frame; b = b->next) {
      this->PrintAlignBestInfo (b);
    }
  }
  else {
    printf ("BESTALIGN ERROR\n");
  }
}


// =======================================================

/*void f() {

  CmdParser cmd(argc, argv);

  cmd.add("-m", "enable XForm and use inXForm      ", "")
     .add("-d", "dir to find hmm definitions       ", "")
     .add("-i", "Output transcriptions to MLF s    ", "")
     .add("-k", "block size for outP calculation   ", "")
     .add("-l", "dir to store label files          ", "")
     .add("-o", "output label formating NCSTWMX    ", "")
     .add("-h", "speaker name pattern              ", "")
     .add("-p", "word insertion penalty            ", "")
     .add("-a", "acoustic scale factor             ", "")
     .add("-r", "pronunciation scale factor        ", "")
     .add("-s", "LM scale factor                   ", "")
     .add("-t", "pruning beam width                ", "")
     .add("-u", "max model pruning                 ", "")
     .add("-v", "wordend beam width                ", "")
     .add("-n", "number of tokens per state        ", "")
     .add("-w", "use language model                ", "")
     .add("-x", "extension for hmm files           ", "")
     .add("-y", "output label file extension       ", "")
     .add("-z", "generate lattices with extension s", "")
     .add("-q", "output lattices format ABtvaldmnr ", "")
     .add("-R", "best align MLF                    ", "")
     .add("-X", "set input lattice extension       ", "")
     .add("-A", "Print command line arguments      ", "")
     .add("-C", "Set config file to cf             ", "")
     .add("-D", "Display configuration variables   ", "")
     .add("-E", "set dir for parent xform to s and optional extension", "")
     .add("-F", "Set source data format to fmt     ", "")
     .add("-H", "Load HMM macro file mmf           ", "")
     .add("-J", "set dir for input xform to s and optional extension", "")
     .add("-L", "Set input label (or net) dir      ", "")
     .add("-P", "Set target label format to fmt    ", "")
     .add("-S", "Set script file to f              ", "")
     .add("-T", "Set trace flags to N              ", "")
     .add("-V", "Print version information         ", "");

  cmd.addGroup(
      "Example: htk/HTKLVRec/HDecode.mod -A -T 1 -a 0.1 -s 1.0 -t 13.0 -z lat"
      "-q tvaldm -o M -l lat -i dev.rec -w data/lm.arpa.txt"
      "-H data/final.mmf -S data/small.scp"
      "data/lexicon.txt data/tiedlist")

  if (!cmd.isOptionLegal())
    cmd.showUsageAndExit();

  // ===========================================================
  hmmDir = cmd["-d"];
  hmmExt = cmd["-x"];
  string mlf_filename = cmd["-i"];
  string target_label_format = cmd["-P"];

  labDir = cmd["-l"];
  labForm = cmd["-o"];
  labExt = cmd["-y"];
  labInExt = cmd["-X"];
  latInDir = cmd["-L"];
  latOutForm = cmd["-q"];
  latOutExt = cmd["-z"];

  insPen = cmd["-p"];

  acScale = cmd["-a"];
  pronScale = cmd["-r"];
  lmScale = cmd["-s"];

  maxModel = cmd["-u"];
  beamWidth = cmd["t"];
  relBeamWidth = beamWidth;

  weBeamWidth = cmd["-v"];
  zsBeamWidth = weBeamWidth;

  langfn = cmd["-w"];

  // ===========================================================
}*/


/*  CC-mode style info for emacs
 Local Variables:
 c-file-style: "htk"
 End:
*/
