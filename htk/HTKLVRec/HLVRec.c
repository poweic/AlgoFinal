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
/*         2002-2003  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HLVRec.h Viterbi recognition engine for       */
/*                        HTK LV Decoder                       */
/* ----------------------------------------------------------- */


char *hlvrec_version = "!HVER!HLVRec:   3.4.1 [GE 12/03/09]";
char *hlvrec_vc_id = "$Id: HLVRec.c,v 1.1.1.1 2006/10/11 09:54:56 jal58 Exp $";


#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HWave.h"
#include "HLabel.h"
#include "HAudio.h"
#include "HParm.h"
#include "HDict.h"
#include "HModel.h"
#include "HUtil.h"
#include "HNet.h"       /* for Lattice -- move to HLattice? */
#include "HAdapt.h"

#include "config.h"

#include "HLVNet.h"
#include "HLVRec.h"
#include "HLVModel.h"

#include <string.h>
#include <assert.h>
#include <misc.h>
#include "HLVRec-propagate.c"
#include "HLVRec-LM.c"

#define PRUNE

/* ----------------------------- Trace Flags ------------------------- */

#define T_TOP 0001         /* Trace  */
#define T_BEST 0002        /* print best token in each frame */
#define T_WORD 0004        /* word end handling */
#define T_TOKSTATS 0010    /* print token/active node stats in each frame */
#define T_PRUNE 0020       /* pruning (need to define DEBUG_TRACE for this!) */
#define T_ACTIV 0040       /* node activation (need to define DEBUG_TRACE for this!) */
#define T_PROP 0100        /* details of token propagation (need to define DEBUG_TRACE for this!) */
#define T_LAT 0200         /* details of lattice generation */
#define T_GC 0400          /* details of garbage collection */
#define T_MEM 01000          /* details of memory usage */

static int trace=0;
static ConfParam *cParm[MAXGLOBS];      /* config parameters */
static int nParm = 0;

static LogFloat maxLMLA = -LZERO; /* maximum jump in LM lookahead per model */

static bool buildLatSE = FALSE;/* build lat from single tok in SENTEND node */

static bool forceLatOut = TRUE;/* always output lattice, even when no token survived */

static int gcFreq = 100;          /* run Garbage Collection every gcFreq frames */

static bool pde = FALSE;      /* partial distance elimination */

static bool useOldPrune = FALSE;     /* backward compatibility for max model and reltok pruning etc. */
static bool mergeTokOnly = TRUE;     /* if merge token set with pruning */
static float maxLNBeamFlr = 0.8;        /* maximum percentile of glogal beam for max model pruning */
static float dynBeamInc = 1.3;          /* dynamic beam increment for max model pruning */
#define LAYER_SIL_NTOK_SCALE 6          /* SIL layer re-adjust token set size e.g. 6 */

/* -------------------------- Global Variables --------------------- */

RelToken startTok = {NULL, NULL, 0.0, 0.0, NULL};
MemHeap recCHeap;                       /* CHEAP for small general allocation */
                                        /* avoid wherever possible! */
static AdaptXForm *inXForm;

/* --------------------------- Prototypes ---------------------- */

/* HLVRec.c */

void InitLVRec(void);
static bool CheckLRTransP (SMatrix transP);
void ReFormatTranscription(Transcription *trans,HTime frameDur,
                           bool states,bool models,bool triStrip,
                           bool normScores,bool killScores,
                           bool centreTimes,bool killTimes,
                           bool killWords,bool killModels);

/* HLVRec-LM.c */
// static void FreeLMCache (LMCache *cache);
// static void ResetLMCache (LMCache *cache);
// LMNodeCache* AllocLMNodeCache (LMCache *cache, int lmlaIdx);

static int LMCacheState_hash (LMState lmstate);
static void PrintPath (DecoderInst *_dec, WordendHyp *we);
static void PrintTok(DecoderInst *_dec, Token *tok);
static void PrintRelTok(DecoderInst *_dec, RelToken *tok);
static void PrintTokSet (DecoderInst *_dec, TokenSet *ts);
TokenSet *BestTokSet (DecoderInst *_dec);
Transcription *TraceBack(MemHeap *heap, DecoderInst *_dec);
static void LatTraceBackCount (DecoderInst *_dec, WordendHyp *path, int *nnodes, int *nlinks);
static void Paths2Lat (DecoderInst *_dec, Lattice *lat, WordendHyp *path,
                       int *na);
Lattice *LatTraceBack (MemHeap *heap, DecoderInst *_dec);
AltWordendHyp *FakeSEpath (DecoderInst *_dec, RelToken *tok, bool useLM);
WordendHyp *AltPathList2Path (DecoderInst *_dec, AltWordendHyp *alt, PronId pron);
WordendHyp *BuildLattice (DecoderInst *_dec);
AltWordendHyp *BuildLatAltList (DecoderInst *_dec, TokenSet *ts, bool useLM);
WordendHyp *BuildForceLat (DecoderInst *_dec);

/* HLVRec-GC.c */
static void MarkModPath (ModendHyp *m);
static void MarkPath (WordendHyp *path);
static void MarkTokSet (TokenSet *ts);
static void SweepPaths (MemHeap *heap);
static void SweepAltPaths (MemHeap *heap);
static void SweepModPaths (MemHeap *heap);

/* HLVRec-outP.c */
static void ResetOutPCache (OutPCache *cache);
static OutPCache *CreateOutPCache (MemHeap *heap, HMMSet *hset, int block);
LogFloat SOutP_ID_mix_Block(HMMSet *hset, int s, Observation *x, StreamElem *se);
void OutPBlock_HMod (StateInfo_lv *si, Observation **obsBlock, 
                     int n, int sIdx, float acScale, LogFloat *outP, int id);

/* HLVRec-misc.c */
void Debug_DumpNet (LexNet *net);

/* --------------------------- Initialisation ---------------------- */


#include "HLVRec-traceback.c"
#include "HLVRec-GC.c"
#include "HLVRec-outP.c"
#include "HLVRec-misc.c"


/* EXPORT->InitLVRec: register module & set configuration parameters */
void InitLVRec(void)
{
   double f;
   int i;
   bool b;
   
   Register (hlvrec_version, hlvrec_vc_id);
   nParm = GetConfig("HLVREC", TRUE, cParm, MAXGLOBS);
   if (nParm>0){
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfFlt (cParm, nParm, "MAXLMLA", &f))
         maxLMLA = f;
      if (GetConfBool (cParm, nParm, "BUILDLATSENTEND",&b)) buildLatSE = b;
      if (GetConfBool (cParm, nParm, "FORCELATOUT",&b)) forceLatOut = b;
      if (GetConfInt (cParm, nParm,"GCFREQ", &i)) gcFreq = i;
      if (GetConfBool (cParm, nParm, "PDE",&b)) pde = b;
      if (GetConfBool (cParm, nParm, "USEOLDPRUNE",&b)) useOldPrune = b;
      if (GetConfBool (cParm, nParm, "MERGETOKONLY",&b)) mergeTokOnly = b;
      if (GetConfFlt (cParm, nParm, "MAXLNBEAMFLR", &f)) maxLNBeamFlr = f;
      if (GetConfFlt (cParm, nParm, "DYNBEAMINC", &f)) dynBeamInc = f;

      if (useOldPrune) {
         mergeTokOnly = FALSE; maxLNBeamFlr = 0.0; dynBeamInc = 1.1;
      }
   }

   CreateHeap (&recCHeap, "Decoder CHEAP", CHEAP, 1, 1.5, 10000, 100000);
}


/* --------------------------- the real code  ---------------------- */


/* CreateDecoderInst

     Create a new instance of the _decInstoding engine. All state information is stored
     here. 
     #### Ideally instances should share other structures (i.e.
          LexNets) this is not implemented, yet.
*/
DecoderInst* Decoder::CreateDecoderInst(HMMSet *hset, FSLM *lm, int nTok, bool latgen, 
                               bool useHModel,
                               int outpBlocksize, bool doPhonePost,
                               bool modAlign)
{
   DecoderInst *_dec;
   int i, N;
   char buf[MAXSTRLEN];

   _dec = (DecoderInst *) New (&recCHeap, sizeof (DecoderInst));

   _dec->lm = lm;
   _dec->hset = hset;
   _dec->useHModel = useHModel;
   /*    _dec->net = net; */

   /* create compact State info. This can change number of shared states! */
   /* #### this is ugly as we end up doing this twice, if we use adaptation! */
   _dec->si = ConvertHSet (&gcheap, hset, _dec->useHModel);

   CreateHeap (&_dec->heap, "Decoder Instance heap", MSTAK, 1, 1.5, 10000, 100000);

   CreateHeap (&_dec->nodeInstanceHeap, "Decoder NodeInstance heap", 
               MHEAP, sizeof (LexNodeInst), 1.5, 1000, 10000);


   _dec->nTok = nTok;
   _dec->latgen = latgen;
   _dec->nLayers = 0;
   _dec->instsLayer = NULL;

   /* alloc & init Heaps for TokenSets */
   N = MaxStatesInSet (_dec->hset);
   _dec->maxNStates = N;

   _dec->tokSetHeap = (MemHeap *) New (&_dec->heap, N * sizeof (MemHeap));


   /* #### make initial size of heap blocks smaller,
      or don't alloc unneeded ones in the first place (scan HMMSet) */
   for (i = 0; i < N; ++i) {
   sprintf (buf, "Decoder %d TokenSet heap", i+1);
      CreateHeap (&_dec->tokSetHeap[i], buf, 
                  MHEAP, (i+1) * sizeof (TokenSet), 9, 10, 5000);
   }   

   _dec->tempTS = (TokenSet **) New (&_dec->heap, N * sizeof (TokenSet *));


   /* alloc Heap for RelToken arrays */
   CreateHeap (&_dec->relTokHeap, "Decoder RelToken array heap",
               MHEAP, _dec->nTok * sizeof (RelToken), 1, 1000, 5000);

   CreateHeap (&_dec->lrelTokHeap, "Decoder RelToken array heap",
               MHEAP, LAYER_SIL_NTOK_SCALE * _dec->nTok * sizeof (RelToken), 1, 1000, 5000);   
   
   /* alloc heap for word end hyps */
   CreateHeap (&_dec->weHypHeap, "WordendHyp heap", MHEAP, sizeof (WordendHyp), 
               1.0, 80000, 800000);
   if (_dec->latgen) {
      CreateHeap (&_dec->altweHypHeap, "AltWordendHyp heap", MHEAP, 
                  sizeof (AltWordendHyp), 1.0, 8000, 80000);
   }
   _dec->modAlign = modAlign;

   if (_dec->modAlign) {
      CreateHeap (&_dec->modendHypHeap, "ModendHyp heap", MHEAP, 
                  sizeof (ModendHyp), 1.0, 80000, 800000);
   }

   /* output probability cache */

   _dec->outPCache = CreateOutPCache (&_dec->heap, _dec->hset, outpBlocksize);

   /* cache debug code */

   /* tag left-to-right models */
   {
      HMMScanState hss;

      NewHMMScan(_dec->hset,&hss);
      do {
         /* #### should check each tidX only once! */
         /*     if (!IsSeenV(hss.hmm->transP)) { */
         if (CheckLRTransP (hss.hmm->transP))
            hss.hmm->tIdx *= -1;
         /*            TouchV(hss.hmm->transP); */
      }
      while(GoNextHMM(&hss));
      EndHMMScan(&hss);
   
   }

   if (doPhonePost)
      InitPhonePost ();
   else
      _dec->nPhone = 0;

   return _dec;
}

/* CheckLRTransP

     determine wheter transition matrix is left-to-right, i.e. no backward transitions
*/
static bool CheckLRTransP (SMatrix transP)
{
   int r,c,N;

   N = NumCols (transP);
   assert (N == NumRows (transP));

   for (r = 1; r <= N; ++r) {
      for (c = 1; c < r; ++c) {
         if (transP[r][c] > LSMALL)
            return FALSE;
      }
      for (c = r+2; c < r; ++c) {
         if (transP[r][c] > LSMALL)
            return FALSE;
      }
   }

   return TRUE;   
}

/* InitDecoderInst

     Initialise previously created _decInstoder instance. This needs to be
     called before each utterance.
*/
void Decoder::InitDecoderInst (LexNet *net, HTime sampRate, LogFloat beamWidth, 
                      LogFloat relBeamWidth, LogFloat weBeamWidth, LogFloat zsBeamWidth,
                      int maxModel, 
                      LogFloat insPen, float acScale, float pronScale, float lmScale,
                      LogFloat fastlmlaBeam)
{       
   int i;

   _dec->net = net;

   if (_dec->nLayers) {
      Dispose (&_dec->heap,_dec->instsLayer);
   }

   /* alloc InstsLayer start pointers */
   _dec->nLayers = net->nLayers;
   _dec->instsLayer = (LexNodeInst **) New (&_dec->heap, net->nLayers * sizeof (LexNodeInst *));

   /* reset inst (i.e. reset pruning, etc.)
      purge all heaps
   */

   ResetHeap (&_dec->nodeInstanceHeap);
   ResetHeap (&_dec->weHypHeap);
   if (_dec->latgen)
      ResetHeap (&_dec->altweHypHeap);
   if (_dec->modAlign)
      ResetHeap (&_dec->modendHypHeap);
   ResetHeap (&_dec->nodeInstanceHeap);
   for (i = 0; i < _dec->maxNStates; ++i) 
      ResetHeap (&_dec->tokSetHeap[i]);
   ResetHeap (&_dec->relTokHeap);
   ResetHeap (&_dec->lrelTokHeap);

   if (trace & T_MEM) {
      printf ("memory stats at start of recognition\n");
      PrintAllHeapStats ();
   }

   _dec->frame = 0;
   _dec->frameDur = sampRate / 1.0e7;
   _dec->maxModel = maxModel;
   _dec->beamWidth = beamWidth;
   _dec->weBeamWidth = weBeamWidth;
   _dec->zsBeamWidth = zsBeamWidth;
   _dec->curBeamWidth = _dec->beamWidth;
   _dec->relBeamWidth = - relBeamWidth;
   _dec->beamLimit = LZERO;

   if (fastlmlaBeam < -LSMALL) {
      _dec->fastlmla = TRUE;
      _dec->fastlmlaBeam = - fastlmlaBeam;
   }
   else {
      _dec->fastlmla = FALSE;
      _dec->fastlmlaBeam = LZERO;
   }

   _dec->tokSetIdCount = 0;

   _dec->insPen = insPen;
   _dec->acScale = acScale;
   _dec->pronScale = pronScale;
   _dec->lmScale = lmScale;

   _dec->maxLMLA = _dec->lmScale * maxLMLA;

   /*      HRec computes interval of possible pre_decInstessor states j for each 
      destination state i in each transition matrix (seIndexes).
   */
   

   /* alloc temp tokenset arrays for use in PropagateInternal */
   for (i=1; i <= _dec->maxNStates; ++i) 
      _dec->tempTS[i] = NewTokSetArrayVar (i, TRUE);

   /* alloc winTok array for MergeTokSet */
   _dec->winTok = (RelToken *) New (&_dec->heap, LAYER_SIL_NTOK_SCALE * _dec->nTok * sizeof (RelToken));


   /* init lists of active LexNode Instances  */
   for (i = 0; i < _dec->nLayers; ++i)
      _dec->instsLayer[i] = NULL;

   /* deactivate all nodes */
   for (i = 0; i < _dec->net->nNodes; ++i) {
      _dec->net->node[i].inst = NULL;
#ifdef COLLECT_STATS_ACTIVATION
      _dec->net->node[i].eventT = -1;
#endif
   }

   ActivateNode (_dec->net->start);
   _dec->net->start->inst->ts[0].n = 1;
   _dec->net->start->inst->ts[0].score = 0.0;
   _dec->net->start->inst->ts[0].relTok[0] = startTok;
   _dec->net->start->inst->ts[0].relTok[0].lmState = LMInitial (_dec->lm);

#ifdef COLLECT_STATS
   _dec->stats.nTokSet = 0;
   _dec->stats.sumTokPerTS = 0;
   _dec->stats.nActive = 0;
   _dec->stats.nActivate = 0;
   _dec->stats.nDeActivate = 0;
   _dec->stats.nFrames = 0;
   _dec->stats.nLMlaCacheHit = 0;
   _dec->stats.nLMlaCacheMiss = 0;
#ifdef COLLECT_STATS_ACTIVATION

   _dec->stats.lnINF = 0;
   {
      int i;
      for (i = 0; i <= STATS_MAXT; ++i)
         _dec->stats.lnDeadT[i] = _dec->stats.lnLiveT[i] = 0;
   }
#endif
#endif

   /* LM lookahead cache */
   _dec->lmCache = CreateLMCache (&_dec->heap);

   /* invalidate OutP cache */
   ResetOutPCache (_dec->outPCache);
}

void Decoder::CleanDecoderInst () {
  _dec->lmCache->free();
}

/* NewTokSetArray

*/
TokenSet *Decoder::NewTokSetArray(int N)
{
   TokenSet *ts;
   int i;

   ts= (TokenSet *) New (&_dec->tokSetHeap[N-1], N * sizeof (TokenSet));

   /* clear token set */
   for (i = 0; i < N; ++i) {
      ts[i].score = 0.0;
      ts[i].n = 0;
      ts[i].id = 0;             /* id=0 means empty TokSet */
      ts[i].relTok = (RelToken *) New (&_dec->relTokHeap, _dec->nTok * sizeof (RelToken));
   }
   return ts;
}

/* NewTokSetArrayVar:

   allocating rel token array for lex node, larger array size
   for silence layer nodes
*/
TokenSet *Decoder::NewTokSetArrayVar(int N, bool isSil)
{
   TokenSet *ts;
   int i;

   ts= (TokenSet *) New (&_dec->tokSetHeap[N-1], N * sizeof (TokenSet));

   /* clear token set */
   for (i = 0; i < N; ++i) {
      ts[i].score = 0.0;
      ts[i].n = 0;
      ts[i].id = 0;             /* id=0 means empty TokSet */
      ts[i].relTok = (isSil) ? (RelToken *) New (&_dec->lrelTokHeap, LAYER_SIL_NTOK_SCALE * _dec->nTok * sizeof (RelToken)) : (RelToken *) New (&_dec->relTokHeap, _dec->nTok * sizeof (RelToken));
   }
   return ts;
}


/* ActivateNode

     Allocate and init new Instance for given node.

*/
LexNodeInst* Decoder::ActivateNode (LexNode *ln)
{
   LexNodeInst *inst;
   int N;               /* number of states in HMM for this node */
   int l;

#ifdef COLLECT_STATS
   ++_dec->stats.nActivate;
#endif
#ifdef COLLECT_STATS_ACTIVATION
   if (ln->eventT != -1) {
      int t;
      t = _dec->frame - ln->eventT;
      if (t > STATS_MAXT)
         t = STATS_MAXT;
      ++_dec->stats.lnDeadT[t];
   }
   ln->eventT = _dec->frame;
#endif

   assert (!ln->inst);

   inst = (LexNodeInst *) New (&_dec->nodeInstanceHeap, 0);

   inst->node = ln;
   ln->inst = inst;

   switch (ln->type) {
   case LN_MODEL:
      N = ln->data.hmm->numStates;
      break;
   case LN_CON:
   case LN_WORDEND:
      N = 1;
      break;
   default:
      abort ();
      break;
   }

   /* alloc N tokensets */
/*    inst->ts = NewTokSetArray (_dec, N);          /\*  size is  N * sizeof (TokenSet) *\/ */

   inst->best = LZERO;

   /* add new instance to list of active nodes in the right place */
   /* find right layer */
   l = _dec->nLayers-1;
   while (_dec->net->layerStart[l] > ln) {
      --l;
      assert (l >= 0);
   }
   inst->ts = NewTokSetArrayVar (N, (l == LAYER_SIL));          /*  size is  N * sizeof (TokenSet) */
#ifdef DEBUG_TRACE
   if (trace & T_ACTIV)
      printf ("allocating %d tokens in array for node in layer %d\n", ((l == LAYER_SIL) ? LAYER_SIL_NTOK_SCALE : 1) * _dec->nTok, l);
#endif

   /* add to linked list */
   inst->next = _dec->instsLayer[l];
   _dec->instsLayer[l] = inst;

#ifdef DEBUG_TRACE
   if (trace & T_ACTIV)
      printf ("activated node in layer %d\n", l);
#endif

   return (inst);
}

void Decoder::DeactivateNode (LexNode *ln)
{
   int N, i, l;

#ifdef COLLECT_STATS
   ++_dec->stats.nDeActivate;
#endif
#ifdef COLLECT_STATS_ACTIVATION
   if (ln->eventT != -1) {
      int t;
      t = _dec->frame - ln->eventT;
      if (t > STATS_MAXT)
         t = STATS_MAXT;
      ++_dec->stats.lnLiveT[t];
   }
   ln->eventT = _dec->frame;
#endif

   assert (ln->inst);
   
   switch (ln->type) {
   case LN_MODEL:
      N = ln->data.hmm->numStates;
      break;
   case LN_CON:
   case LN_WORDEND:
      N = 1;
      break;
   default:
      abort ();
      break;
   }
   
#if 1
   /* find right layer */
   l = _dec->nLayers-1;
   while (_dec->net->layerStart[l] > ln) {
      --l;
      assert (l >= 0);
   }
   for (i = 0; i < N; ++i) {
      if (l == LAYER_SIL) Dispose (&_dec->lrelTokHeap, ln->inst->ts[i].relTok);
      else Dispose (&_dec->relTokHeap, ln->inst->ts[i].relTok);
   }

   Dispose (&_dec->tokSetHeap[N-1], ln->inst->ts);
   Dispose (&_dec->nodeInstanceHeap, ln->inst);
#endif

   ln->inst = NULL;
}


/* PruneTokSet

     apply global and relative beams to a tokenset

     #### this is rather inefficient
*/
void Decoder::PruneTokSet (TokenSet *ts)
{
   RelTokScore deltaLimit;
   RelToken *tok, *dest;
   int i, newN;

   return;

   /* main and relative beam pruning */
   deltaLimit = std::max(_dec->beamLimit - ts->score, _dec->relBeamWidth);
   
   if (deltaLimit > 0) {        /* prune complete TokeSet */
      ts->n = 0;
      ts->id = 0;
      return;
   }

   /* #### maybe don't perform relTok pruning to keep relTokID the same?  */

   newN = 0;            /* number of relToks kept */
   for (i = 0, dest = tok = ts->relTok; i < ts->n; ++i, ++tok) {
      if (tok->delta > deltaLimit) {   /* keep */
         //if (dest != tok)
	 *dest = *tok;
         ++dest;
         ++newN;
      }
   }

   /* #### could calculate newN from difference between dest and ts->tok */
   if (newN != ts->n) {         /* some RelToks got pruned! */
      ts->n = newN;
      ts->id = ++_dec->tokSetIdCount;

   }
}


/* stolen from HRec.c */

/* EXPORT->FormatTranscription: Format transcription prior to output */
void ReFormatTranscription(Transcription *trans,HTime frameDur,
                         bool states,bool models,bool triStrip,
                         bool normScores,bool killScores,
                         bool centreTimes,bool killTimes,
                         bool killWords,bool killModels)
{
   LabList *ll;
   LLink lab;
   HTime end;
   char buf[MAXSTRLEN],*p,tail[64];
   int lev,j,frames;
   
   if (killScores) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
            lab->score=0.0;
            for (j=1;j<=ll->maxAuxLab;j++)
               lab->auxScore[j]=0.0;
         }
      }
   }
   if (triStrip) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
            if (states && !models) {
               strcpy(buf,lab->labid->name);
               if ((p=strrchr(buf,'['))!=NULL) {
                  strcpy(tail,p);
                  *p=0;
               }
               else
                  *tail=0;
               TriStrip(buf); strcat(buf,tail);
               lab->labid=GetLabId(buf,TRUE);
            }
            else {
               strcpy(buf,lab->labid->name);
               TriStrip(buf); lab->labid=GetLabId(buf,TRUE);
            }
            for (j=1;j<=ll->maxAuxLab;j++) {
               if (lab->auxLab[j]==NULL) continue;
               strcpy(buf,lab->auxLab[j]->name);
               TriStrip(buf); lab->auxLab[j]=GetLabId(buf,TRUE);
            }
         }
      }
   }
   if (normScores) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
            frames=(int)floor((lab->end-lab->start)/frameDur + 0.4);
            if (frames==0) lab->score=0.0;
            else lab->score=lab->score/frames;
            if (states && models && ll->maxAuxLab>0 && lab->auxLab[1]!=NULL) {
               end=AuxLabEndTime(lab,1);
               frames=(int)floor((end-lab->start)/frameDur + 0.4);
               if (frames==0) lab->auxScore[1]=0.0;
               else lab->auxScore[1]=lab->auxScore[1]/frames;
            }
         }
      }
   }
   if (killTimes) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
            lab->start=lab->end=-1.0;
         }
      }
   }
   if (centreTimes) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
            lab->start+=frameDur/2;
            lab->end-=frameDur/2;
         }
      }
   }
   if (killWords) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         if (ll->maxAuxLab>0)
            for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ)
               lab->auxLab[ll->maxAuxLab]=NULL;
      }
   }
   if (killModels && models && states) {
      for (lev=1;lev<=trans->numLists;lev++) {
         ll=GetLabelList(trans,lev);
         if (ll->maxAuxLab==2)
            for(lab=ll->head->succ;lab->succ!=NULL;lab=lab->succ) {
               lab->auxLab[1]=lab->auxLab[2];
               lab->auxScore[1]=lab->auxScore[2];
               lab->auxLab[2]=NULL;
            }
      }
   }
}

/*  CC-mode style info for emacs
 Local Variables:
 c-file-style: "htk"
 End:
*/
