#include <misc.h>

/* 
 * external global variables 
 * */
extern int trace;

/* 
 * class HiddenMarkovModel
 * */

HiddenMarkovModel::HiddenMarkovModel() {
  /* init model heap & set early to support loading MMFs */
  CreateHeap(&_modelHeap, "Model heap",  MSTAK, 1, 0.0, 100000, 800000 );
  CreateHMMSet(&_hset, &_modelHeap, true); 
}

void HiddenMarkovModel::init(char* hmmList_fn, char* mmf_fn, char* inXFormDir_fn,
    char* hmmDir, char* hmmExt) {

  if (MakeHMMSet (&_hset, hmmList_fn) < SUCCESS) 
    HError (4128, "Initialise: MakeHMMSet failed");

  AddMMF (&_hset, mmf_fn); 

  if (inXFormDir_fn)
    AddInXFormDir(&_hset, inXFormDir_fn);

  /* Read accoustic models */
  if (trace & T_TOP)
    printf ("Reading acoustic models...");

  // hmmDir and hmmExt seems to be NULL
  if (LoadHMMSet (&_hset, hmmDir, hmmExt) < SUCCESS) 
    HError (4128, "Initialise: LoadHMMSet failed");

  /* convert to INVDIAGC */
  ::ConvDiagC (&_hset, true);
  ::ConvLogWt (&_hset);

  if (trace&T_TOP)
    printf("Read %d physical / %d logical HMMs\n", _hset.numPhyHMM, _hset.numLogHMM);  
}

void HiddenMarkovModel::SetStreamWidths(bool *eSep) {
  ::SetStreamWidths(_hset.pkind, _hset.vecSize, _hset.swidth, eSep);
}

bool HiddenMarkovModel::UpdateSpkrStats(XFInfo *xfinfo, char *datafn) {
  return ::UpdateSpkrStats(&_hset, xfinfo, datafn);
}

MLink HiddenMarkovModel::FindMacroStruct(char type, Ptr structure) {
  return ::FindMacroStruct (&_hset, type, structure);
}

MLink HiddenMarkovModel::FindMacroName(char type, LabId id) {
  return ::FindMacroName(&_hset, type, id);
}

HMMSet& HiddenMarkovModel::get_hset() { return _hset; }

/*
 * class Vocabulary
 * */

Vocabulary::Vocabulary():
  startWord("<s>"), endWord("</s>"), spModel("sp"), silModel("sil"), silDict(false) {
}

void Vocabulary::init(char* dictFn) {
  dict_fn = dictFn;

  /* Read dictionary */
  if (trace & T_TOP)
    printf ("Reading dictionary from %s\n", dict_fn);

  InitVocab (&_vocab);
  if (ReadDict (dict_fn, &_vocab) < SUCCESS)
    HError (9999, "Initialise: ReadDict failed");
}

void Vocabulary::Process() {
   /* process dictionary */
   startLab = GetLabId (startWord, false);
   if (!startLab) 
      HError (9999, "HDecode: cannot find STARTWORD '%s'\n", startWord);

   endLab = GetLabId (endWord, false);
   if (!endLab) 
      HError (9999, "HDecode: cannot find ENDWORD '%s'\n", endWord);

   spLab = GetLabId (spModel, false);
   if (!spLab)
      HError (9999, "HDecode: cannot find label 'sp'");

   silLab = GetLabId (silModel, false);
   if (!silLab)
      HError (9999, "HDecode: cannot find label 'sil'");

}

void Vocabulary::ConvertSilDict () {
  _vocab.ConvertSilDict(spLab, silLab, startLab, endLab);
}

void Vocabulary::MarkAllProns () {
  _vocab.MarkAllProns();
}

void Vocabulary::MarkAllWords () {
  _vocab.MarkAllWords();
}

void Vocabulary::UnMarkAllWords () {
  _vocab.UnMarkAllWords();
}

void Vocabulary::MarkAllWordsfromLat (Lattice* lat) {
  ::MarkAllWordsfromLat(&_vocab, lat, silDict);
}

void Vocabulary::CheckForSkipInSpModel() {

}

Vocab& Vocabulary::get_vocab() {
  return _vocab;
}

/*
 * Class LanguageModel 
 * */

LanguageModel::LanguageModel(Vocabulary& vocab): lm(NULL), _vocab(vocab) {
  CreateHeap (&_lmHeap, "LM heap", MSTAK, 1, 0,1000000, 10000000);
}

void LanguageModel::loadFromFile(char* lm_fn) {

  if (trace & T_TOP)
    printf ("Reading language model from %s\n", lm_fn);
      
  if (!lm_fn) 
    HError (9999, "HDecode: no LM or lattice specified");

  lm = CreateLM (&_lmHeap, lm_fn, _vocab.startWord, _vocab.endWord, &_vocab.get_vocab());
}

void LanguageModel::loadFromLattice(char* lat_fn, Lattice* lat) {
  lm = CreateLMfromLat (&_lmHeap, lat_fn, lat, &_vocab.get_vocab());
}

FSLM* LanguageModel::get_lm() {
  return lm;
}

Lattice* LanguageModel::ReadLattice(FILE *file, bool shortArc, bool add2Dict) {
  return ::ReadLattice(file, &_lmHeap, &_vocab.get_vocab(), shortArc, add2Dict);
}

void LanguageModel::ResetHeap() {
  ::ResetHeap(&_lmHeap);
}

Decoder::Decoder(LanguageModel& lm, Vocabulary& vocab, HiddenMarkovModel& hmm, XFInfo& xfInfo)
  : net(NULL), _lm(lm), _vocab(vocab), _hmm(hmm), _dec(NULL), _xfInfo(xfInfo),

    /* default configs */
    insPen(0.0), acScale(1.0), pronScale(1.0), lmScale(1.0), maxModel(0),
    beamWidth(-LZERO), weBeamWidth(-LZERO), zsBeamWidth(-LZERO), relBeamWidth(-LZERO),
    latPruneBeam(-LZERO), latPruneAPS(0), fastlmlaBeam(-LZERO),
    nTok(32), useHModel(false), outpBlocksize(1), obs(NULL), 

    bestAlignMLF(NULL),

    /* other configs */
    ofmt(UNDEFF), labDir(NULL), labExt("rec"), labForm(NULL), 
    latRescore(false), latInDir(NULL), latInExt("lat"), latFileMask(NULL), 
    latGen(false), latOutDir(NULL), latOutExt("lat"), latOutForm(NULL), 
    dataForm(UNDEFF), mergeTokOnly(true), gcFreq(100), dynBeamInc(1.3) {

      // Nothing to do now.
}

void Decoder::init() {

   /* init Heaps */
   CreateHeap (&_transHeap,"Transcription heap",MSTAK,1,0,8000,80000);

   _vocab.Process();

   if (_vocab.silDict) {    /* dict contains -/sp/sil variants (with probs) */
      _vocab.ConvertSilDict();

      /* check for skip in sp model */
      LabId spLab = GetLabId ("sp", false);
      if (!spLab)
	HError (9999, "cannot find 'sp' model.");

      MLink spML = _hmm.FindMacroName('l', spLab);
      if (!spML)
	HError (9999, "cannot find model for sp");

      HLink spHMM = (HLink) spML->structure;
      int N = spHMM->numStates;

      if (spHMM->transP[1][N] > LSMALL)
	HError (9999, "HDecode: using -/sp/sil dictionary but sp contains tee transition!");
   }
   else {       /* lvx-style dict (no sp/sil at wordend */
     _vocab.MarkAllProns();
   }

   if (!latRescore) {

      /* mark all words  for inclusion in Net */
     _vocab.MarkAllWords();

      /* create network */
      net = new LexNet;
      net->init(&_vocab.get_vocab(), &_hmm.get_hset(), _vocab.startWord, _vocab.endWord, _vocab.silDict);
      
      /* Read language model */
      _lm.loadFromFile(_lm.langfn);
   }

   bool modAlign = false;
   if (latOutForm) {
      if (strchr (latOutForm, 'd'))
         modAlign = true;
      if (strchr (latOutForm, 'n'))
         HError (9999, "DoRecognition: likelihoods for model alignment not supported");
   }

   /* create Decoder instance */
   _dec = CreateDecoderInst (&_hmm.get_hset(), _lm.get_lm(), nTok, true, useHModel, outpBlocksize,
                            bestAlignMLF ? true : false,
                            modAlign);
   
   /* create buffers for observations */
   bool eSep;
   _hmm.SetStreamWidths(&eSep);

   obs = (Observation *) New (&gcheap, outpBlocksize * sizeof (Observation));
   for (int i = 0; i < outpBlocksize; ++i)
      obs[i] = MakeObservation (&gcheap, _hmm.get_hset().swidth, _hmm.get_hset().pkind, 
                                (_hmm.get_hset().hsKind == DISCRETEHS), eSep);

   CreateHeap (&_inputBufHeap, "Input Buffer Heap", MSTAK, 1, 1.0, 80000, 800000);

   /* Initialise adaptation */

   /* sort out masks just in case using adaptation */
   if (_xfInfo.inSpkrPat == NULL) _xfInfo.inSpkrPat = _xfInfo.outSpkrPat; 
   if (_xfInfo.paSpkrPat == NULL) _xfInfo.paSpkrPat = _xfInfo.outSpkrPat; 

   if (_xfInfo.useOutXForm) {
      CreateHeap(&_regHeap,   "regClassStore",  MSTAK, 1, 0.5, 1000, 8000 );
      /* This initialises things - temporary hack - THINK!! */
      CreateAdaptXForm(&_hmm.get_hset(), "tmp");

      /* online adaptation not supported yet! */
   }
}

void Decoder::rescoreLattice(char* fn) {
  /* read lattice and create LM */
  char latfn[MAXSTRLEN],buf3[MAXSTRLEN];
  FILE *latF;
  bool isPipe;
  Lattice *lat;

  /* clear out previous LexNet, Lattice and LM structures */
  _lm.ResetHeap();

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
    lat = _lm.ReadLattice(latF, false, false);
    // lat = ReadLattice (latF, &_lmHeap, &_vocab.get_vocab(), false, false);
    FClose (latF, isPipe);
    if (!lat)
      HError (9999, "DoRecognition: cannot read lattice file %s\n", latfn);
  }

  /* mark prons of all words in lattice */
  _vocab.UnMarkAllWords();
  _vocab.MarkAllWordsfromLat(lat);

  /* create network of all the words/prons marked (word->aux and pron->aux == 1) */
  if (trace & T_TOP)
    printf ("Creating network\n");
  net->init(&_vocab.get_vocab(), &_hmm.get_hset(), _vocab.startWord, _vocab.endWord, _vocab.silDict);
  // net = CreateLexNet (&_netHeap, &_vocab.get_vocab(), &_hmm.get_hset(), _vocab.startWord, _vocab.endWord, _vocab.silDict);

  /* create LM based on pronIds defined by CreateLexNet */
  if (trace & T_TOP)
    printf ("Creating language model\n");

  _lm.loadFromLattice(latfn, lat);
  _dec->lm = _lm.get_lm();
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
      bool changed;
      changed = _hmm.UpdateSpkrStats(&_xfInfo, fn);
   }

   startClock = clock();

   /* get transcrition of 1-best alignment */
   if (bestAlignMLF)
      bestAlignInfo = CreateBestInfo (fn, pbInfo.tgtSampRate/1.0e7);
   
   parmBuf = OpenBuffer (&_inputBufHeap, fn, 50, dataForm, TRI_UNDEF, TRI_UNDEF);
   if (!parmBuf)
      HError (9999, "HDecode: Opening input failed");
   
   GetBufferInfo (parmBuf, &pbInfo);
   if (pbInfo.tgtPK != _hmm.get_hset().pkind)
      HError (9999, "HDecode: Incompatible parm kinds %s vs. %s",
              ParmKind2Str (pbInfo.tgtPK, buf1),
              ParmKind2Str (_hmm.get_hset().pkind, buf2));
              
   if (latRescore)
     rescoreLattice(fn);

   if (weBeamWidth > beamWidth)
      weBeamWidth = beamWidth;
   if (zsBeamWidth > beamWidth)
      zsBeamWidth = beamWidth;

   InitDecoderInst (net, pbInfo.tgtSampRate, beamWidth, relBeamWidth,
                    weBeamWidth, zsBeamWidth, maxModel,
                    insPen, acScale, pronScale, lmScale, fastlmlaBeam);

   net->vocabFN = _vocab.dict_fn;
   _dec->utterFN = fn;

   frameN = frameProc = 0;
   while (BufferStatus (parmBuf) != PB_CLEARED) {
      ReadAsBuffer (parmBuf, &obs[frameN % outpBlocksize]);

      if (frameN+1 >= outpBlocksize) {  /* enough frames available */
         if (trace & T_OBS)
            PrintObservation (frameProc+1, &obs[frameProc % outpBlocksize], 13);
         for (i = 0; i < outpBlocksize; ++i)
            obsBlock[i] = &obs[(frameProc + i) % outpBlocksize];

#ifdef DEBUG_TRACE
         // fprintf(stdout, "\nProcessing frame %d :\n", frameProc);
#endif
         
         ProcessFrame (obsBlock, outpBlocksize, _xfInfo.inXForm);
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
      
      ProcessFrame (obsBlock, bs, _xfInfo.inXForm);
      if (bestAlignInfo)
         this->AnalyseSearchSpace (bestAlignInfo);
      ++frameProc;
   }
   assert (frameProc == frameN);

   
   endClock = clock();
   cpuSec = (endClock - startClock) / (double) CLOCKS_PER_SEC;
   printf ("CPU time %f  utterance length %f  RT factor %f\n",
           cpuSec, frameN*_dec->frameDur, cpuSec / (frameN*_dec->frameDur));

   trans = TraceBack (&_transHeap);

   /* save 1-best transcription */
   /* the following is from HVite.c */
   if (trans) {
      char labfn[MAXSTRLEN];

      if (labForm != NULL)
         ReFormatTranscription (trans, pbInfo.tgtSampRate, false, false,
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
      lat = LatTraceBack (&_transHeap);

      /* prune lattice */
      if (lat && latPruneBeam < - LSMALL) {
         lat = LatPrune (&_transHeap, lat, latPruneBeam, latPruneAPS);
      }

      /* the following is from HVite.c */
      if (lat) {
         char latfn[MAXSTRLEN];
         char *p;
         bool isPipe;
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


   printf ("Stats: nTokSet %lu\n", _dec->stats.nTokSet);
   printf ("Stats: TokPerSet %f\n", _dec->stats.sumTokPerTS / (double) _dec->stats.nTokSet);
   printf ("Stats: activePerFrame %f\n", _dec->stats.nActive / (double) _dec->stats.nFrames);
   printf ("Stats: activateNodePerFrame %f\n", _dec->stats.nActivate / (double) _dec->stats.nFrames);
   printf ("Stats: deActivateNodePerFrame %f\n\n", 
           _dec->stats.nDeActivate / (double) _dec->stats.nFrames);
#ifdef COLLECT_STATS_ACTIVATION
   {
      int i;
      for (i = 0; i <= STATS_MAXT; ++i)
         printf ("T %d Dead %lu Live %lu\n", i, _dec->stats.lnDeadT[i], _dec->stats.lnLiveT[i]);
   }
#endif

   if (trace & T_MEM) {
      printf ("memory stats at end of recognition\n");
      PrintAllHeapStats ();
   }

   ResetHeap (&_inputBufHeap);
   ResetHeap (&_transHeap);
   CleanDecoderInst ();
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
   if (bestTrans->head->tail->pred->pred->labid == _vocab.spLab ||
       bestTrans->head->tail->pred->pred->labid == _vocab.silLab) {
      LLink delLL;
      
      delLL = bestTrans->head->tail->pred->pred;
      /* add sp's frames (if any) to final sil */
      delLL->succ->start = delLL->pred->end;
      
      delLL->pred->succ = delLL->succ;
      delLL->succ->pred = delLL->pred;
   }
   
   ln = net->start;
   assert (ln->type == LN_MODEL);
   m = FindMacroStruct (&_hmm.get_hset(), 'h', ln->data.hmm);
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

/**********  align best code  ****************************************/

/* find the LN_MODEL lexnode following ln that has label lab
   step over LN_CON and LN_WORDEND nodes.
   return NULL if not found
*/
BestInfo* Decoder::FindLexNetLab (MemHeap *heap, LexNode *ln, LLink ll, HTime frameDur) {
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
      phonePost = _dec->phonePost[(int) monoPhone->aux];
   } else
      phonePost = 999.99;

   int l = _dec->nLayers-1;
   while (_dec->net->layerStart[l] > b->ln) {
      --l;
      assert (l >= 0);
   }
   
   printf ("BESTALIGN frame %4d best %.3f alignbest %d -> %d ln %p layer %d score %.3f phonePost %.3f\n", 
           _dec->frame, _dec->bestScore, 
           b->start, b->end, b->ln, l, score, phonePost);
}

void Decoder::AnalyseSearchSpace (BestInfo *bestInfo) {
  if (!bestInfo)
    return;

  BestInfo *b;
  LabId monoPhone;

  monoPhone =(LabId) _dec->bestInst->node->data.hmm->hook;
  printf ("frame %4d best %.3f phonePost %.3f\n", _dec->frame, 
      _dec->bestScore, _dec->phonePost[(int) monoPhone->aux]);

  for (b = bestInfo; b; b = b->next) {
    if (b->start < _dec->frame && b->end >= _dec->frame) 
      break;
  }
  if (b) {
    this->PrintAlignBestInfo (b);
    for (b = b->next; b && b->start == b->end && b->start == _dec->frame; b = b->next) {
      this->PrintAlignBestInfo (b);
    }
  }
  else {
    printf ("BESTALIGN ERROR\n");
  }
}
