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

    void SetStreamWidths(bool *eSep);

    bool UpdateSpkrStats(XFInfo *xfinfo, char *datafn);

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

    bool silDict;
};

class LanguageModel {
  public:
    LanguageModel(Vocabulary& vocab);

    void loadFromFile(char* lm_fn);

    void loadFromLattice(char* lat_fn, Lattice* lat);

    FSLM* get_lm();

    Lattice *ReadLattice(FILE *file, bool shortArc, bool add2Dict);

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

    //==== These functions are moved here by coldsheep 0613

    // From propagate.c
  
    void HandleWordend (LexNode *ln);
    void PropagateExternal (LexNodeInst *inst, bool handleWE, bool wintTree);
    void UpdateWordEndHyp (LexNodeInst *inst);
    void AddPronProbs (TokenSet *ts, int var);
    void MergeTokSet (TokenSet *src, TokenSet *dest, LogFloat score, bool prune);

    void TokSetBucketSortPruning(TokenSet *dest, const RelTokScore& deltaLimit,
	int nWinTok, RelToken* &winTok, TokScore &winScore);

    void FindWinningTokens(TokenSet *src, TokenSet *dest, LogFloat score,
	RelTokScore srcCorr, RelTokScore destCorr, RelTokScore deltaLimit,
	RelToken* &winTok, int &nWinTok, int* nWin);

    void PropagateInternal ();
    void PropagateInternal (LexNodeInst *inst);
    void PropIntoNode (TokenSet *ts, LexNode *ln, bool updateLMLA);
    void UpdateModPaths (TokenSet *ts, LexNode *ln);
    void HandleSpSkipLayer (LexNodeInst *inst);
    void ProcessFrame (Observation **obsBlock, int nObs, AdaptXForm *xform);

    void SetObservation (Observation **obsBlock, int nObs);
    void WordEndBeamPruning (LexNodeInst* head, TokScore &beamLimit);
    void ZSLayerBeamPruning (LexNodeInst* head, TokScore &beamLimit);
    void RelaxBeamLimit ();

    void __collect_stats__ (TokenSet *instTS, int N); 
    void OptLeftToRightPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm);
    // void GeneralPropagateInternal();
    void GeneralPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm);

    // From HLVRec-LM.c
    LMTokScore LMLA_nocache (LMState lmState, int lmlaIdx);
    LMCache* CreateLMCache (MemHeap *heap);
    LMTokScore LMCacheLookaheadProb (LMState lmState, int lmlaIdx, bool fastlmla);
    LMTokScore LMCacheTransProb (FSLM *lm, LMState src, PronId pronid, LMState *dest);
    void UpdateLMlookahead(LexNode *ln);

    // From HLVRec.c 3/7 are moved here
    DecoderInst *CreateDecoderInst(HMMSet *hset, FSLM *lm, int nTok, bool latgen, 
	  bool useHModel, int outpBlocksize, bool doPhonePost, bool modAlign);

    void InitDecoderInst ( LexNet *net, HTime sampRate, LogFloat beamWidth, 
	  LogFloat relBeamWidth, LogFloat weBeamWidth, LogFloat zsBeamWidth,
	  int maxModel, 
	  LogFloat insPen, float acScale, float pronScale, float lmScale,
	  LogFloat fastlmlaBeam);

    void CleanDecoderInst ();
    TokenSet *NewTokSetArray(int N);
    TokenSet *NewTokSetArrayVar(int N, bool isSil);

    LexNodeInst *ActivateNode (LexNode *ln);
    void DeactivateNode (LexNode *ln);
    void PruneTokSet (TokenSet *ts);

    // From outP.c
    LogFloat cOutP (Observation *x, HLink hmm, int state);

    // From GC.c
    void GarbageCollectPaths ();

    // Form HLVRec-misc.c
    void CheckTokenSetOrder (TokenSet *ts);
    void CheckTokenSetId (TokenSet *ts1, TokenSet *ts2);
    WordendHyp *CombinePaths (RelToken *winner, RelToken *loser, LogFloat diff);
    void Debug_Check_Score ();
    void InitPhonePost ();
    void CalcPhonePost ();
    void AccumulateStats ();

    // From traceback.c 13/18 are moved here
    AltWordendHyp* FakeSEpath (RelToken *tok, bool useLM);
    AltWordendHyp* BuildLatAltList (TokenSet *ts, bool useLM);
    WordendHyp *BuildForceLat ();
    WordendHyp* BuildLattice ();
    Lattice* LatTraceBack (MemHeap *heap);
    void PrintPath (WordendHyp *we);
    void PrintTok(Token *tok);
    void PrintRelTok(RelToken *tok);
    void PrintTokSet (TokenSet *ts);
    TokenSet* BestTokSet ();
    Transcription* TraceBack(MemHeap *heap);
    void LatTraceBackCount (WordendHyp *path, int *nnodes, int *nlinks);
    void Paths2Lat (Lattice *lat, WordendHyp *path, int *na);
    LAlign* LAlignFromModpath (MemHeap *heap, ModendHyp *modpath, int wordStart, short *nLAlign);
    LAlign *LAlignFromAltModpath (MemHeap *heap, ModendHyp *modpath, ModendHyp *mainModpath, int wordStart, short *nLAlign);
    void PrintModPath (ModendHyp *m);
    void CheckLAlign (Lattice *lat);

    // ===== THIS function does NOT belongs here. Move it somewhere else. =====
    void rescoreLattice(char* fn);

  private:
    LexNet* net;
    LanguageModel& _lm;
    Vocabulary& _vocab;
    HiddenMarkovModel& _hmm;
    XFInfo& _xfInfo;

    DecoderInst* _decInst;

    MemHeap _inputBufHeap;
    MemHeap _transHeap;
    MemHeap _regHeap;
  public:

    LogFloat insPen;   /* word insertion penalty */

    float acScale;     /* acoustic scaling factor */
    float pronScale;   /* pronunciation scaling factor */
    float lmScale;     /* LM scaling factor */

    int maxModel;	    /* max model pruning */
    LogFloat beamWidth;     /* pruning global beam width */
    LogFloat weBeamWidth;   /* pruning wordend beam width */
    LogFloat zsBeamWidth;   /* pruning z-s beam width */
    LogFloat relBeamWidth;  /* pruning relative beam width */
    LogFloat latPruneBeam;  /* lattice pruning beam width */
    LogFloat latPruneAPS;   /* lattice pruning arcs per sec limit */
    LogFloat fastlmlaBeam;  /* do fast LM la outside this beam */

    int nTok;		    /* number of different LMStates per HMM state */
    bool useHModel;	    /* use standard HModel OutP functions */
    int outpBlocksize;	    /* number of frames for which outP is calculated in one go */
    Observation *obs;       /* array of Observations */

    /* info for comparing scores from alignment of 1-best with search */
    char *bestAlignMLF;     /* MLF with 1-best alignment */

    FileFormat ofmt;	    /* Label output file format */
    char *labDir;	    /* output label file directory */
    char *labExt;	    /* output label file extension */
    char *labForm;	    /* output label format */

    bool latRescore;	    /* read lattice for each utterance and rescore? */
    char *latInDir;	    /* lattice input directory */
    char *latInExt;	    /* latttice input extension */
    char *latFileMask;	    /* mask for reading lattice */

    bool latGen;	    /* output lattice? */
    char *latOutDir;	    /* lattice output directory */
    char *latOutExt;	    /* latttice output extension */
    char *latOutForm;	    /* lattice output format */

    FileFormat dataForm;    /* data input file format */
    bool mergeTokOnly;	    /* if merge token set with pruning */
    AdaptXForm *inXForm;
    float dynBeamInc;       /* dynamic beam increment for max model pruning */
    int gcFreq;             /* run Garbage Collection every gcFreq frames */
    Stats stats;            /* statistics about pruning etc. */
};

