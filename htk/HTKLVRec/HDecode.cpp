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


#include <cmdparser.h>
#include "misc.h"

int trace = 0;

/* -------------------------- Global Variables etc ---------------------- */

// char *langfn;		/* LM filename from commandline */
char *dictfn;		/* dict filename from commandline */
char *hmmListfn;		/* model list filename from commandline */
char *hmmDir = NULL;     /* directory to look for HMM def files */
char *hmmExt = NULL;     /* HMM def file extension */

LexNet *net = NULL;      /* Lexicon network of all required words/prons */

XFInfo xfInfo;

/* -------------------------- Prototypes -------------------------------- */

void ParseCommandArguments(Vocabulary& vocab, LanguageModel& lm, HiddenMarkovModel& hmm, Decoder& decoder);
void InitAll(int argc, char *argv[]);

void SetConfParms (Vocabulary& vocab, Decoder& decoder);
void ReportUsage ();
Boolean UpdateSpkrModels (char *fn);

/* ---------------- Configuration Parameters ---------------------------- */

static ConfParam *cParm[MAXGLOBS];
static int nParm = 0;		/* total num params */


/* ---------------- Debug support  ------------------------------------- */

#if 0
FILE *debug_stdout = stdout;
FILE *debug_stderr = stderr;
#endif

/* ---------------- Process Command Line ------------------------- */

char *mmf_fn = NULL;
char *inXFormDir_fn = NULL;

int main (int argc, char *argv[]) {

  // ===== Create Vocabulary, LM, Lexicon, HMM, Decoder =====
  Vocabulary vocab;
  LanguageModel lm(vocab);
  HiddenMarkovModel hmm;
  Decoder decoder(lm, vocab, hmm, xfInfo);

  // ===== Settings =====
  InitAll(argc, argv);
  ParseCommandArguments(vocab, lm, hmm, decoder);

  vocab.init(dictfn);
  hmm.init(hmmListfn, mmf_fn, inXFormDir_fn, hmmDir, hmmExt);
  decoder.init();

  SetConfParms(vocab, decoder);

  /* load 1-best alignment */
  if (decoder.bestAlignMLF)
    LoadMasterFile (decoder.bestAlignMLF);

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
  hmm.UpdateSpkrStats(&xfInfo, NULL);

  Exit(0);             /* maybe print config and exit */
  return 0;
}

/* SetConfParms: set conf parms relevant to this tool */
void SetConfParms (Vocabulary& vocab, Decoder& decoder) {
   int i;
   double f;
   Boolean b;
   char buf[MAXSTRLEN];

   nParm = GetConfig ("HDECODE", TRUE, cParm, MAXGLOBS);
   if (nParm > 0) {
      if (GetConfInt (cParm, nParm, "TRACE", &i))	  trace			= i;
      if (GetConfStr (cParm, nParm, "STARTWORD", buf))	  vocab.startWord	= CopyString (&gstack, buf);
      if (GetConfStr (cParm, nParm, "ENDWORD", buf))	  vocab.endWord		= CopyString (&gstack, buf);
      if (GetConfFlt (cParm, nParm, "LATPRUNEBEAM", &f))  decoder.latPruneBeam  = f;
      if (GetConfFlt (cParm, nParm, "FASTLMLABEAM", &f))  decoder.fastlmlaBeam  = f;
      if (GetConfFlt (cParm, nParm, "LATPRUNEAPS", &f))	  decoder.latPruneAPS	= f;
      if (GetConfStr (cParm, nParm, "BESTALIGNMLF", buf)) decoder.bestAlignMLF	= CopyString (&gstack, buf);
      if (GetConfBool (cParm, nParm, "USEHMODEL",&b))	  decoder.useHModel	= b;
      if (GetConfStr(cParm,nParm,"LATFILEMASK",buf))	  decoder.latFileMask	= CopyString(&gstack, buf);
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


/**********  align best code  ****************************************/

/* find the LN_MODEL lexnode following ln that has label lab
   step over LN_CON and LN_WORDEND nodes.
   return NULL if not found
*/
BestInfo* Decoder::FindLexNetLab (MemHeap *heap, LexNode *ln, LLink ll, HTime frameDur)
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
      m = _hmm.FindMacroStruct('h', follLN->data.hmm);
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
}

void ParseCommandArguments(Vocabulary& vocab, LanguageModel& lm, HiddenMarkovModel& hmm, Decoder& decoder) {
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
	 if ((decoder.ofmt = Str2Format (GetStrArg ())) == ALIEN)
	    HError (-3289,
		    "HVite: Warning ALIEN Label output file format set");
	 break;

      case 'l':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Label/Lattice output directory expected");
	 decoder.labDir = GetStrArg ();
         decoder.latOutDir = decoder.labDir;
	 break;
      case 'o':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output label format expected");
	 decoder.labForm = GetStrArg ();
	 break;
      case 'y':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output label file extension expected");
	 decoder.labExt = GetStrArg ();
	 break;

      case 'X':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Input Lattice file extension expected");
	 decoder.latInExt = GetStrArg ();
	 break;
      case 'L':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Input Lattice directory expected");
	 decoder.latInDir = GetStrArg ();
	 break;

      case 'q':
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output lattice format expected");
	 decoder.latOutForm = GetStrArg ();
	 break;
      case 'z':
         decoder.latGen = TRUE;
	 if (NextArg () != STRINGARG)
	    HError (4019, "HDecode: Output lattice file extension expected");
	 decoder.latOutExt = GetStrArg ();
	 break;

      case 'p':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: word insertion penalty expected");
         decoder.insPen = GetFltArg ();
         if (decoder.insPen > 0.0)
            HError (-1, "HDecode: positive word insertion penalty???");
	 break;
      case 'a':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: acoustic scale factor expected");
	 decoder.acScale = GetFltArg ();
	 break;
      case 'r':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: pronunciation scale factor expected");
	  decoder.pronScale = GetFltArg ();
          vocab.silDict = TRUE;       /* #### maybe separate switch for this? */
	 break;
      case 's':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: LM scale factor expected");
         decoder.lmScale= GetFltArg ();
	 break;


      case 'u':
	 if (NextArg () != INTARG)
	    HError (4019, "HDecode: max model pruning limit expected");
         decoder.maxModel = GetIntArg ();
	 break;

      case 't':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: beam width expected");
	 decoder.beamWidth = GetFltArg ();
         if (decoder.latPruneBeam == -LZERO)
            decoder.latPruneBeam = decoder.beamWidth;
         decoder.relBeamWidth = decoder.beamWidth;
         if (NextArg () == FLOATARG)
            decoder.relBeamWidth = GetFltArg ();
	 break;

      case 'v':
	 if (NextArg () != FLOATARG)
	    HError (4019, "HDecode: wordend beam width expected");
         decoder.weBeamWidth = GetFltArg ();
         decoder.zsBeamWidth = decoder.weBeamWidth;
	 if (NextArg () == FLOATARG)
            decoder.zsBeamWidth = GetFltArg ();
         break;

      case 'w':
	 if (NextArg() != STRINGARG) {
            /*	    HError (4119,"HDecode: LM File name expected"); */
            decoder.latRescore = TRUE;
         }
         else
            lm.langfn = GetStrArg();
	 break;

      case 'n':
	 decoder.nTok = GetChkedInt (0, 1024, s);
	 break;

      case 'k':
	 decoder.outpBlocksize = GetChkedInt (0, MAXBLOCKOBS, s);
	 break;

      case 'H':
	 if (NextArg() != STRINGARG)
	    HError (4119,"HDecode: MMF File name expected");
	 
	 mmf_fn = GetStrArg();
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
	 inXFormDir_fn = GetStrArg();
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
	 decoder.bestAlignMLF = GetStrArg ();
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

   if (decoder.beamWidth > -LSMALL)
      HError (4019, "main beam is too wide!");

   if (xfInfo.useInXForm) {
      if (!decoder.useHModel) {
         HError (-4019, "HDecode: setting USEHMODEL to TRUE.");
         decoder.useHModel = TRUE;
      }
      if (decoder.outpBlocksize != 1) {
         HError (-4019, "HDecode: outP blocksize >1 not supported with new XForm code! setting to 1.");
         decoder.outpBlocksize = 1;
      }
   }   
}

/*****************  main recognition function  ************************/

Boolean UpdateSpkrModels (char *fn) {
   HError (1, "MLLR or FV transforms not supported");
   return FALSE;
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
