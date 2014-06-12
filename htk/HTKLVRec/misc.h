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

/* -------------------------- Trace Flags & Vars ------------------------ */

#define T_TOP 00001		/* Basic progress reporting */
#define T_OBS 00002		/* Print Observation */
#define T_ADP 00004		/* Adaptation */
#define T_MEM 00010		/* Memory usage, start and finish */

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

class HiddenMarkovModel {
  public:
    HiddenMarkovModel();

    void init(char* hmmList_fn, char* mmf_fn, char* inXFormDir_fn,
	char* hmmDir, char* hmmExt);

    void SetStreamWidths(Boolean *eSep);

    Boolean UpdateSpkrStats(XFInfo *xfinfo, char *datafn);

    MLink FindMacroStruct(char type, Ptr structure);

    MLink FindMacroName(char type, LabId id);

    HMMSet& get_hset();

  private:
    MemHeap _modelHeap;
    HMMSet _hset;
};

class Vocabulary {
  public:
    Vocabulary();

    void init(char* dict_fn);

    void Process();

    void ConvertSilDict ();

    void MarkAllProns ();

    void MarkAllWords ();

    void UnMarkAllWords ();

    void MarkAllWordsfromLat (Lattice* lat);

    void CheckForSkipInSpModel();

    Vocab& get_vocab();

    void set_start_word(char* str) {
      startWord = str;
    }

    void set_end_word(char* str) {
      endWord = str;
    }

  private:
    Vocab _vocab;

  public:
    char* dict_fn;

    char *startWord; /* word used at start of network */
    LabId startLab;  /*   corresponding LabId */
    char *endWord;   /* word used at end of network */
    LabId endLab;    /*   corresponding LabId */

    char *spModel;   /* model used as word end Short Pause */
    LabId spLab;     /*   corresponding LabId */
    char *silModel;  /* model used as word end Silence */
    LabId silLab;    /*   corresponding LabId */

    Boolean silDict;
};

class LanguageModel {
  public:
    LanguageModel(Vocabulary& vocab);

    void loadFromFile(char* lm_fn);

    void loadFromLattice(char* lat_fn, Lattice* lat);

    FSLM* get_lm();

    Lattice *ReadLattice(FILE *file, Boolean shortArc, Boolean add2Dict);

    void ResetHeap();

  private:
    FSLM *lm;         /* language model */
    MemHeap _lmHeap;

    Vocabulary& _vocab;

  public:
    char *langfn;     /* LM filename from commandline */
};

class Decoder {
  public:
    Decoder(LanguageModel& lm, Vocabulary& vocab, HiddenMarkovModel& hmm, XFInfo& xfInfo);
    void init();
    void recognize(char* fn);

    BestInfo* CreateBestInfo (char* fn, HTime frameDur);
    BestInfo* FindLexNetLab (MemHeap *heap, LexNode *ln, LLink ll, HTime frameDur);

    void PrintAlignBestInfo (BestInfo *bestInfo);
    void AnalyseSearchSpace (BestInfo *bestInfo);


    // ===== THIS function does NOT belongs here. Move it somewhere else. =====
    void rescoreLattice(char* fn);

  private:
    LanguageModel& _lm;
    Vocabulary& _vocab;
    HiddenMarkovModel& _hmm;
    XFInfo& _xfInfo;

    DecoderInst* _decInst;

    MemHeap _netHeap;
    // MemHeap _lmHeap;
    MemHeap _inputBufHeap;
    MemHeap _transHeap;
    MemHeap _regHeap;
  public:

    LogFloat insPen;   /* word insertion penalty */

    float acScale;     /* acoustic scaling factor */
    float pronScale;   /* pronunciation scaling factor */
    float lmScale;     /* LM scaling factor */

    int maxModel;        /* max model pruning */
    LogFloat beamWidth;     /* pruning global beam width */
    LogFloat weBeamWidth;   /* pruning wordend beam width */
    LogFloat zsBeamWidth;   /* pruning z-s beam width */
    LogFloat relBeamWidth;  /* pruning relative beam width */
    LogFloat latPruneBeam;  /* lattice pruning beam width */
    LogFloat latPruneAPS;   /* lattice pruning arcs per sec limit */
    LogFloat fastlmlaBeam;  /* do fast LM la outside this beam */

    int nTok;		    /* number of different LMStates per HMM state */
    Boolean useHModel;	    /* use standard HModel OutP functions */
    int outpBlocksize;	    /* number of frames for which outP is calculated in one go */
    Observation *obs;       /* array of Observations */

    /* info for comparing scores from alignment of 1-best with search */
    char *bestAlignMLF;     /* MLF with 1-best alignment */

    FileFormat ofmt;	    /* Label output file format */
    char *labDir;	    /* output label file directory */
    char *labExt;	    /* output label file extension */
    char *labForm;	    /* output label format */

    Boolean latRescore;     /* read lattice for each utterance and rescore? */
    char *latInDir;	    /* lattice input directory */
    char *latInExt;	    /* latttice input extension */
    char *latFileMask;	    /* mask for reading lattice */

    Boolean latGen;	    /* output lattice? */
    char *latOutDir;	    /* lattice output directory */
    char *latOutExt;	    /* latttice output extension */
    char *latOutForm;	    /* lattice output format */

    FileFormat dataForm;    /* data input file format */
};

