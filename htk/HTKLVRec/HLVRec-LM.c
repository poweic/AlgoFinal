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
/*         File: HLVRec-LM.c Update LM lookahead for           */
/*                           HTK LV Decoder                    */
/* ----------------------------------------------------------- */

/* UpdateLMlookahead

     update LM lookahead info for all RelToks in entry state of ln

     This might change the best score in ts-score, but will not affect the
     tokenset identity ts->id. All tokensets we might later want to merge with
     will be in the same node and have thus undergone the same LMla change.
*/
void Decoder::UpdateLMlookahead(LexNode *ln)
{
  LMTokScore lmscore;

  TokenSet *ts = ln->inst->ts;
  RelToken *tok	= ts->relTok;
  RelTokScore bestDelta = LZERO;
  unsigned int lmlaIdx = ln->lmlaIdx;

  for (int i = 0; i < ts->n; ++i, ++tok) {

    if (!_decInst->fastlmla)
      lmscore = LMCacheLookaheadProb (tok->lmState, lmlaIdx, FALSE);
    else {    /* if we ever do fastLMLA, be careful as tok->lmscore might increase! */
      lmscore = LMCacheLookaheadProb (tok->lmState, lmlaIdx, tok->delta < _decInst->fastlmlaBeam);
      lmscore = std::min(lmscore, tok->lmscore);
    }

    if (lmscore > LSMALL &&  tok->lmscore - lmscore > _decInst->maxLMLA)
      lmscore = tok->lmscore - _decInst->maxLMLA;

    tok->delta += lmscore - tok->lmscore;     /* subtract previous lookahead */
    bestDelta = std::max(bestDelta, tok->delta);

    tok->lmscore = lmscore;   /* store current lookahead */
  }

  /* renormalise to new best score */

  if (bestDelta > LSMALL) {
    tok = ts->relTok;
    for (int i = 0; i < ts->n; ++i, ++tok)
      tok->delta -= bestDelta;
    ts->score += bestDelta;
  }
  else {       /* short cut pruning for LMLA = LZERO */
    ts->n = 0;
    ts->score = LZERO;
    ts->id = 0;
  }
}


/******************* LM trans * lookahead caching */

LMCache* Decoder::CreateLMCache (MemHeap *heap)
{
   LMCache *cache;
   int i;

   cache = (LMCache *) New (heap, sizeof (LMCache));
   CreateHeap (&cache->nodeHeap, "LMNodeCache Heap", MHEAP, sizeof (LMNodeCache),
               1.0, 1000, 2000);

   cache->nNode = _decInst->net->laTree->nNodes + _decInst->net->laTree->nCompNodes;
   cache->node = (LMNodeCache **) New (heap, cache->nNode * sizeof (LMNodeCache *));
   for (i = 0; i < cache->nNode; ++i)
      cache->node[i] = NULL;

   cache->transHit = cache->transMiss = 0;
   cache->laHit = cache->laMiss = 0;
   return cache;
}

static int LMCacheState_hash (LMState lmstate)
{
   return ((unsigned int) lmstate % LMCACHE_NLA);
}

/* LMCacheTransProb

     return the (scaled!) LM transition prob and dest LMState for the given LMState and PronId
*/
LMTokScore Decoder::LMCacheTransProb (FSLM *lm, LMState src, PronId pronid, LMState *dest)
{
   return _decInst->lmScale * LMTransProb (lm, src, pronid, dest);
}

LMTokScore Decoder::LMLA_nocache (LMState lmState, int lmlaIdx)
{
   LMlaTree *laTree;
   LMTokScore lmscore;

   laTree = _decInst->net->laTree;
   if (lmlaIdx < laTree->nNodes) {        /* simple node */
      LMlaNode *laNode;
         
      laNode = &laTree->node[lmlaIdx];
      lmscore = _decInst->lmScale * LMLookAhead (_decInst->lm, lmState, 
                                            laNode->loWE, laNode->hiWE);
   }
   else {         /* complex node */
      CompLMlaNode *laNode;
      LMTokScore score;
      int i;
         
      laNode = &laTree->compNode[lmlaIdx - laTree->nNodes];
      
      lmscore = LZERO;
      for (i = 0; i < laNode->n; ++i) {
         score = LMLA_nocache (lmState, laNode->lmlaIdx[i]);
         if (score > lmscore)
            lmscore = score;
      }
   }
   return lmscore;
}

/* LMCacheLookaheadProb

     return the (scaled!) LM lookahead score for the given LMState and lmlaIdx
*/
LMTokScore Decoder::LMCacheLookaheadProb (LMState lmState, int lmlaIdx, bool fastlmla)
{
   LMCache *cache;
   LMTokScore lmscore;
   LMNodeCache *nodeCache;
   LMCacheLA *entry;
   int i;

   cache = _decInst->lmCache;
   assert (lmlaIdx < cache->nNode);
   nodeCache = cache->node[lmlaIdx];

   if (fastlmla) {      /* #### should only go to fast LMState if real one is not cahced */
      lmState = Fast_LMLA_LMState (_decInst->lm, lmState);
   }

   if (nodeCache) {  
      /* touch this state */
      nodeCache->t = _decInst->frame;
      
      for (i = 0; i < nodeCache->nEntries; ++i) {
         entry = &nodeCache->la[i];
         if (entry->src == lmState)
            break;
      }
      if (i < nodeCache->nEntries) {
         ++cache->laHit;
         return entry->prob;
      }
      entry = &nodeCache->la[nodeCache->nextFree];
      nodeCache->nextFree = (nodeCache->nextFree + 1) % nodeCache->size;
      if (nodeCache->nEntries < nodeCache->size)
         ++nodeCache->nEntries;
   }
   else {
      /* alloc new LMNodeCache */
      nodeCache = cache->node[lmlaIdx] = cache->alloc(lmlaIdx);

      entry = &nodeCache->la[0];
   }
   
   /* now entry points to the place to store the prob we are about to calulate */
   {
      LMlaTree *laTree;
      
      laTree = _decInst->net->laTree;
      if (lmlaIdx < laTree->nNodes) {        /* simple node */
         LMlaNode *laNode;
         
         laNode = &laTree->node[lmlaIdx];
         ++cache->laMiss;
         lmscore = _decInst->lmScale * LMLookAhead (_decInst->lm, lmState, 
                                               laNode->loWE, laNode->hiWE);
      }
      else {         /* complex node */
         CompLMlaNode *laNode;
         LMTokScore score;
         int i;
         
         laNode = &laTree->compNode[lmlaIdx - laTree->nNodes];
         
         lmscore = LZERO;
         for (i = 0; i < laNode->n; ++i) {
            score = LMCacheLookaheadProb (lmState, laNode->lmlaIdx[i], FALSE);
            if (score > lmscore)
               lmscore = score;
         }
      }
   }

   if (lmscore < LSMALL)
      lmscore = LZERO;
   
   entry->src = (void**) lmState;
   entry->prob = lmscore;

   /*    printf ("lmla %f\n", lmscore); */
   return lmscore;
}

