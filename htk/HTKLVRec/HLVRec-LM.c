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
   int i;
   TokenSet *ts;
   RelToken *tok;
   unsigned int lmlaIdx;
   LMTokScore lmscore;
   RelTokScore bestDelta;

   lmlaIdx = ln->lmlaIdx;
   ts = ln->inst->ts;

   bestDelta = LZERO;

   for (i = 0, tok = ts->relTok; i < ts->n; ++i, ++tok) {

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
      for (i = 0, tok = ts->relTok; i < ts->n; ++i, ++tok)
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

static void FreeLMCache (LMCache *cache)
{
   DeleteHeap (&cache->nodeHeap);
}

static void ResetLMCache (LMCache *cache)
{
   int i;
   
   ResetHeap (&cache->nodeHeap);
   for (i = 0; i < cache->nNode; ++i)
      cache->node[i] = NULL;

   cache->transHit = cache->transMiss = 0;
   cache->laHit = cache->laMiss = 0;
}

static int LMCacheState_hash (LMState lmstate)
{
   return ((unsigned int) lmstate % LMCACHE_NLA);
}

#if 0
static int LMCacheTrans_hash (PronId pron)
{
   return ((unsigned int) pron % LMCACHE_NTRANS);
}

static int LMCacheLA_hash (int idx)
{
   return ((unsigned int) idx % LMCACHE_NLA);
}
#endif

LMNodeCache* AllocLMNodeCache (LMCache *cache, int lmlaIdx)
{
   LMNodeCache *n;

#if 0
   printf ("new LMNodeCache %d\n", lmlaIdx);
#endif
   n = (LMNodeCache *) New (&cache->nodeHeap, sizeof (LMNodeCache));
   memset ((void *) n, 1, sizeof (LMNodeCache));  /* clear all src entries */
   n->idx = lmlaIdx;
   n->size = LMCACHE_NLA;
   n->nextFree = n->nEntries = 0;

   return n;
}

/* LMCacheTransProb

     return the (scaled!) LM transition prob and dest LMState for the given LMState and PronId
*/
LMTokScore Decoder::LMCacheTransProb (FSLM *lm, LMState src, PronId pronid, LMState *dest)
{
   return _decInst->lmScale * LMTransProb (lm, src, pronid, dest);

#if 0
   LMCache *cache;
   int hash;
   LMStateCache *stateCache;
   LMCacheTrans *entry;

   cache = _decInst->lmCache;
   hash = LMStateCache_hash (src);

   for (stateCache = cache->state[hash]; stateCache; stateCache = stateCache->next)
      if (stateCache->src == src)
         break;
   
   if (stateCache) {  
      /* touch this state */
      stateCache->t = _decInst->frame;
      
      entry = &stateCache->trans[LMCacheTrans_hash (pronid)];
      
      if (entry->pronid == pronid) {
         ++cache->transHit;
         *dest = entry->dest;
         return entry->prob;
      }
   }
   else {
      /* alloc new LMStateCache */
      stateCache = AllocLMStateCache (cache, src);
      stateCache->next = cache->state[hash];
      cache->state[hash] = stateCache;

      entry = &stateCache->trans[LMCacheTrans_hash (pronid)];
   }

   /* now entry points to the place to store the prob we are about to calulate */

   ++cache->transMiss;
   entry->pronid = pronid;
   entry->prob = _decInst->lmScale * LMTransProb (lm, src, pronid, dest);
   entry->dest = *dest;
   return entry->prob;
#endif
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
#if 0         /* #### very expensive sanity check */
         assert (entry->prob == LMLA_nocache (_decInst, lmState, lmlaIdx));
#endif
         return entry->prob;
      }
      entry = &nodeCache->la[nodeCache->nextFree];
      nodeCache->nextFree = (nodeCache->nextFree + 1) % nodeCache->size;
      if (nodeCache->nEntries < nodeCache->size)
         ++nodeCache->nEntries;
   }
   else {
      /* alloc new LMNodeCache */
      nodeCache = cache->node[lmlaIdx] = AllocLMNodeCache (cache, lmlaIdx);

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

