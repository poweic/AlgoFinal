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
/*         File: HLVRec-propagate.c Viterbi recognition engine */
/*                                  for HTK LV Decoder, token  */
/*                                  propagation                */
/* ----------------------------------------------------------- */

#include <cmath>
#include <algorithm>

char *hlvrec_prop_vc_id = "$Id: HLVRec-propagate.c,v 1.1.1.1 2006/10/11 09:54:56 jal58 Exp $";

/* perform Bucket sort/Histogram pruning to reduce to _decInst->nTok tokens */
void Decoder::TokSetBucketSortPruning(TokenSet *dest, const RelTokScore& deltaLimit,
    int nWinTok, RelToken* &winTok, TokScore &winScore) {

  const size_t NBINS = 64;

  int i, j;

  dest->id = ++_decInst->tokSetIdCount;    /* #### new id always necessary? */
  LogFloat binWidth = deltaLimit * 1.001 / NBINS;   /* handle delta == deltaLimit case */

  int n[NBINS] = {0};

  for (i = 0; i < nWinTok; ++i)
    ++n[(int) (winTok[i].delta / binWidth)];

  int nTok = 0;
  i = -1;
  while (nTok < _decInst->nTok)
    nTok += n[++i];

  if (nTok == _decInst->nTok) {
    int binLimit = i;

    for (i = 0, j = 0; i < nWinTok; ++i) { 
      if ((int) (winTok[i].delta / binWidth) <= binLimit) {
	dest->relTok[j] = winTok[i];
	++j;
      }
    }
  }
  else {
    int nBetter;
    LogFloat bestDelta;

    /* do not include last bin */
    LogFloat limit = binWidth * i;
    nTok -= n[i]; 

    /* need to relax limit so that we get an extra (_decInst->nTok - nTok) tokens */
    /* #### very simplistic implementation -- imporve? */

    bestDelta = limit;
    do {
      limit = bestDelta;
      bestDelta = LZERO;
      nBetter = 0;
      for (i = 0, j = 0; i < nWinTok; ++i) {
	if (winTok[i].delta >= limit)
	  ++nBetter;
	else if (winTok[i].delta > bestDelta)
	  bestDelta = winTok[i].delta;
      }
    } while (nBetter < _decInst->nTok);

    /* multiple tokens with delta == limit ==> delete some */
    if (nBetter > _decInst->nTok) {  
      for (i = 0; nBetter > _decInst->nTok; ++i)
	if (winTok[i].delta == limit) {
	  winTok[i].delta = LZERO;
	  --nBetter;
	}
    }

    for (i = 0, j = 0; i < nWinTok; ++i) { 
      if (winTok[i].delta >= limit) {
	dest->relTok[j] = winTok[i];
	++j;
      }
    }
  }

  dest->n = _decInst->nTok;
  dest->score = winScore;

}

void Decoder::FindWinningTokens(TokenSet *src, TokenSet *dest, LogFloat score,
    RelTokScore srcCorr, RelTokScore destCorr, RelTokScore deltaLimit,
    RelToken* &winTok, int &nWinTok, int* nWin) {

  RelToken *srcTok = src->relTok, *destTok = dest->relTok;
  int nSrcTok = src->n, nDestTok = dest->n;
  int winLoc;
  do {
    if (TOK_LMSTATE_EQ(srcTok, destTok)) {
      /* pick winner */
      if (src->score + srcTok->delta + score > dest->score + destTok->delta) {
	/* store srcTok */
	winTok[nWinTok] = *srcTok;
	winTok[nWinTok].delta += srcCorr + score;
	winLoc = 0;
      } else {
	/* store destTok */
	winTok[nWinTok] = *destTok;
	winTok[nWinTok].delta += destCorr;
	winLoc = 1;
      }
      ++srcTok;
      --nSrcTok;
      ++destTok;
      --nDestTok;
    } else if (TOK_LMSTATE_LT(srcTok, destTok)) {
      /* store srcTok */
      winTok[nWinTok] = *srcTok;
      winTok[nWinTok].delta += srcCorr + score;
      winLoc = 0;
      ++srcTok;
      --nSrcTok;
    } else {
      /* store destTok */
      winTok[nWinTok] = *destTok;
      winTok[nWinTok].delta += destCorr;
      winLoc = 1;
      ++destTok;
      --nDestTok;
    }

    if (winTok[nWinTok].delta >= deltaLimit) {      /* keep or prune? */
      ++nWinTok;
      ++nWin[winLoc];
    }

  } while (nSrcTok != 0 && nDestTok != 0);

  /* Add left overs to winTok set, only at most one of the two loops will
   * actually do something */
  for (int i = nSrcTok; i > 0; --i, ++srcTok) {
    winTok[nWinTok] = *srcTok;
    winTok[nWinTok].delta += srcCorr + score;
    if (winTok[nWinTok].delta >= deltaLimit) {     /* keep or prune? */
      ++nWinTok;
      ++nWin[0];
    }
  }

  for (int i = nDestTok; i > 0; --i, ++destTok) {
    winTok[nWinTok] = *destTok;
    winTok[nWinTok].delta += destCorr;
    if (winTok[nWinTok].delta >= deltaLimit) {     /* keep or prune? */
      ++nWinTok;
      ++nWin[1];
    }
  }
}

/* MergeTokSet

     Merge TokenSet src into dest after adding score to all src scores
     by recombining tokens in the same LM state and keeping only the
     _decInst->nTok best tokens in different LM states.

*/
void Decoder::MergeTokSet (TokenSet *src, TokenSet *dest, LogFloat score, bool prune) {
   /* empty dest tokenset -> just copy src */
   if (dest->n == 0) {   
      *dest = *src;
      dest->score += score;

   } else if (prune && src->score + score < _decInst->beamLimit) {
     // Do nothing
   }
#ifdef TSIDOPT
   else if (src->id == dest->id) {      /* TokenSet Id optimisation from [Odell:2000] */
      /* Only compare Tokensets' best scores and pick better */
      /* id, n and all RelToks are the same anyway */
      dest->score = std::max(src->score + score, dest->score);
   }
#endif
   else {    
     /* expensive MergeTokSet, #### move into separate function */
     /* exploit & retain RelTok order (sorted by lmState?) */
     TokScore winScore;
     RelTokScore srcCorr, destCorr, deltaLimit;

     /* first go at sorted Tok merge, very explicit, no optimisation at all, yet! */

     /* find best score */
     if (src->score + score > dest->score) {
       winScore = src->score + score;
       srcCorr = - score;     /* avoid adding score twice! */
       destCorr = dest->score - winScore;
     } else {
       winScore = dest->score;
       srcCorr = src->score - winScore;
       destCorr = 0.0;
     }

     /* scaled relative beam, must initialize !!!*/;
     deltaLimit = _decInst->nTok * _decInst->relBeamWidth;  

     /* set pruning deltaLimit */
     if (prune)
       deltaLimit = std::max(_decInst->beamLimit - winScore, _decInst->relBeamWidth);

     /* find winning tokens */
     /* location where winnning toks came from: 0 == src, 1 == dest */
     RelToken *winTok = _decInst->winTok;
     int nWinTok = 0;
     int nWin[2] = {0, 0};
     FindWinningTokens(src, dest, score, srcCorr, destCorr, deltaLimit, winTok, nWinTok, nWin);

     if (nWinTok <= _decInst->nTok) {
       /* just copy */
       for (int i = 0; i < nWinTok; ++i)
	 dest->relTok[i] = winTok[i];

       dest->n = nWinTok;
       dest->score = winScore;

       if (nWin[0] == nWinTok)
	 dest->id = src->id;          /* copy src->id */
       else if (nWin[1] == nWinTok)
	 dest->id = dest->id;         /* copy dest->id */
       else
	 dest->id = ++_decInst->tokSetIdCount;    /* new id */

     } else 
       TokSetBucketSortPruning(dest, deltaLimit, nWinTok, winTok, winScore);
   }
}

void Decoder::__collect_stats__ (TokenSet *instTS, int N) {
#ifdef COLLECT_STATS
  for (int i = 0; i < N; ++i) {
    if (instTS[i].n > 0) {
      ++_decInst->stats.nTokSet;
      _decInst->stats.sumTokPerTS += instTS[i].n;
    }
  }
#endif
}

static int PI_LR = 0;
static int PI_GEN = 0;

void Decoder::OptLeftToRightPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm) {

  PI_LR++;
  LogFloat bestScore = LZERO;

  /* loop transition for state N-1 (which has no forward trans) */
  instTS[N-2].score += trP[N-1][N-1];

  for (int i = N-2; i >=2; --i) {             /* for all internal states, except N-1 */
    /* ts is an array of tokenSet ( of the type TokenSet* ) ,
     * starting from state 0 to stae N ( N = 5 typically ). */
    TokenSet* ts = &instTS[i-1];

    if (ts->n > 0) {
      /* propagate forward from state i */

      /* only propagate to next -- no skip!! ( i.e. i->i+1 no i->i+2 ) */
      MergeTokSet (ts, ts+1 , trP[i][i+1], (!mergeTokOnly));

      /* loop transition */
      ts->score += trP[i][i];
    }
  }

  /* entry transition i=1 -> j=2 */
  if ((instTS[0].n > 0) && (trP[1][2] > LSMALL)) {
    MergeTokSet (&instTS[0], &instTS[1], trP[1][2], (!mergeTokOnly));
  }

  /* output probabilities */
  for (int i = 2; i < N; ++i) {             /* for all internal states */
    if (instTS[i-1].n > 0) {
      LogFloat outP = cOutP (_decInst->obs, hmm, i);
      instTS[i-1].score += outP;         /* only change top score */
      if (instTS[i-1].score > bestScore)
	bestScore = instTS[i-1].score;
    }
  }
  inst->best = bestScore;

  instTS[0].n = 0;             /* clear entry state */
  instTS[0].id = 0;

  instTS[N-1].n = 0;           /* clear exit state */
  instTS[N-1].id = 0;

  /* exit transition i=N-1 -> j=N */
  /* don't prune -- beam stil refers to last frame's scores, here we have 
     already added this frame's outP */
  if (instTS[N-2].n > 0)
    MergeTokSet (&instTS[N-2], &instTS[N-1], trP[N-1][N], FALSE);

  if (bestScore > _decInst->bestScore) {
    _decInst->bestScore = bestScore;
    _decInst->bestInst = inst;
  }

  __collect_stats__(instTS, N);
}

void Decoder::GeneralPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm) {

  /* temp storage for N tokensets  */
  TokenSet *tempTS = _decInst->tempTS[N];

  PI_GEN++;

  /* internal propagation; transition i -> j,  \forall 2 <= j <= N-1 */

  /* internal states */
  for (int j = 2; j < N; ++j) {
    tempTS[j-1].score = 0.0;
    tempTS[j-1].n = 0;
    tempTS[j-1].id = 0;

    for (int i = 1; i < N; ++i) {
      if ((instTS[i-1].n > 0) && (trP[i][j] > LSMALL))
	MergeTokSet (&instTS[i-1], &tempTS[j-1], trP[i][j], (!mergeTokOnly));
    }

    if (tempTS[j-1].n > 0) {
      LogFloat outP = cOutP (_decInst->obs, hmm, j);
      tempTS[j-1].score += outP;         /* only change top tok */
    }
  }

  instTS[0].n = 0;          /* clear entry state */
  instTS[0].id = 0;

  /* internal states: copy temp array back and find best score */
  LogFloat bestScore = LZERO;
  for (int j = 1; j < N-1; ++j) {
    instTS[j] = tempTS[j]; // copy reltoks

    if (instTS[j].n > 0)
      bestScore = std::max(instTS[j].score, bestScore);
  }
  inst->best = bestScore;

  /* exit state (j=N), merge directly into ints->ts */
  int j = N;
  instTS[j-1].n = 0;
  instTS[j-1].id = 0;

  for (int i = 2; i < N; ++i) {
    if ((instTS[i-1].n > 0) && (trP[i][j] > LSMALL))
      MergeTokSet (&instTS[i-1], &instTS[j-1], trP[i][j], FALSE);
  }

  /* exit state score is some internal state score plus transP,
     thus can be ignored for findeing the best score */

  /* update global best score */
  if (bestScore > _decInst->bestScore) {
    _decInst->bestScore = bestScore;
    _decInst->bestInst = inst;
  }

  /* # this only collects stats for the model nodes */
  __collect_stats__(instTS, N);
}

/* PropagateInternal
     Internal token propagation
*/
void Decoder::PropagateInternal (LexNodeInst *inst) {

   LexNode *ln = inst->node;
   HLink hmm = ln->data.hmm;
   int N = hmm->numStates;
   SMatrix trP = hmm->transP;
   TokenSet* instTS = inst->ts;

   assert (ln->type == LN_MODEL);               /* Model node */
   /* LM lookahead has already been updated in PropIntoNode() !!! */

   /* Main beam pruning: prune tokensets before propagation
      the beamLimit is the one found during the last frame */
   auto ts = &instTS[0];
   for (int i = 1; i < N; ++i, ++ts) {
     if (ts->score < _decInst->beamLimit)
       ts->n = ts->id = 0;
   }

   /* Optimised version for L-R models OR 
      General (not left-to-right) Propagate Internal   */
   if (hmm->tIdx < 0)
      OptLeftToRightPropagateInternal(inst, instTS, N, trP, hmm);
   else
     GeneralPropagateInternal(inst, instTS, N, trP, hmm);
}

#ifdef MODALIGN
void Decoder::UpdateModPaths (TokenSet *ts, LexNode *ln) {
   if (ln->type == LN_CON)
      return;

   ModendHyp *m;
   RelToken *tok;
   int i;
   
   /* don't accumulate info for CON nodes */
   /* #### optimise by sharing ModendHyp's between tokens with
      same tok->modpath */
   for (i = 0, tok = ts->relTok; i < ts->n; ++i, ++tok) {
      m = (ModendHyp*) New (&_decInst->modendHypHeap, sizeof (ModendHyp));
      m->frame = _decInst->frame;
      m->ln = ln;
      m->prev = tok->modpath;

      tok->modpath = m;
   }
}
#endif   

/* PropIntoNode
     Propagate tokenset into entry state of LexNode, activating as necessary
*/
void Decoder::PropIntoNode (TokenSet *ts, LexNode *ln, bool updateLMLA)
{
   LexNodeInst *inst;
   TokScore best;

   if (!ln->inst)                /* activate if necessary */
      ActivateNode (ln);

   inst = ln->inst;
         
   /* propagate tokens from ln's exit into follLN's entry state */
   MergeTokSet (ts, &inst->ts[0], 0.0, TRUE);
   
   /* only update inst->best if no LMLA update necessary or already done! */
   if (ln->type == LN_WORDEND || ln->lmlaIdx == 0){
      inst->best = std::max(inst->ts[0].score, inst->best);
      return;
   }

   if (updateLMLA) {
      UpdateLMlookahead (ln);
      inst->best = std::max(inst->ts[0].score, inst->best);
   }
}

/* PropagateExternal

     External token propagation. Activate following nodes if necessary and propagate 
     token set into their entry states.
*/
void Decoder::PropagateExternal ( LexNodeInst *inst, bool handleWE, bool wintTree)
{
   LexNode *ln, *follLN;
   int i, N;
   TokenSet *entryTS, *exitTS;
   
   ln = inst->node;
   entryTS = &inst->ts[0];

   /* handle tee transition and set N */
   if (ln->type == LN_MODEL) {
      HLink hmm;

      hmm = ln->data.hmm;
      N = hmm->numStates;
      exitTS = &inst->ts[N-1];

      /* main beam pruning: only propagate if above beamLimit
         #### maybe unnecessary? can be done in PropIntoNode */
      if ((entryTS->n > 0) && (hmm->transP[1][N] > LSMALL) && 
          (entryTS->score > _decInst->beamLimit)) {
         MergeTokSet (entryTS, exitTS, hmm->transP[1][N], TRUE);
      }
   }
   else {
      N = 1;    /* all other nodes are 1 state */
      exitTS = entryTS;

      /* for wordend nodes first apply LM and recombine */
      if (ln->type == LN_WORDEND) {
         /* for LAYER_WR call HandleWordend() on all WE nodes in separate pass for WE-pruning */
         if (handleWE)
            HandleWordend (ln);
      }      
      /* main beam pruning: only propagate if above beamLimit */
      if (exitTS->score < _decInst->beamLimit) {
         exitTS->n = 0;
         exitTS->id = 0;
         inst->best = LZERO;
      }
      else
         inst->best = exitTS->score;
   }

   /* prune token sets in states 2..N-1 */
   /* ####RELTOK  don't do this to avoid changing tokSet->id?? */
   for (i = 2; i < N; ++i) {
      if (inst->ts[i-1].n > 0)
         PruneTokSet (&inst->ts[i-1]);
   }
   /* prune exit state token set */
   // the following two lines are not used in our L0 L1 test.
   if (exitTS->n > 0)
      PruneTokSet (exitTS);


   /* any tokens in exit state? */
   if (exitTS->n > 0 && exitTS->score > _decInst->beamLimit) {
#ifdef MODALIGN
      if (_decInst->modAlign)
         UpdateModPaths (exitTS, ln);
#endif
      /* loop over following nodes */
      for (i = 0; i < ln -> nfoll; ++i) {
         follLN = ln->foll[i];
         PropIntoNode (exitTS, follLN, wintTree);
      } /* for i following nodes */
   }
}


/* HandleWordend
     update traceback, add LM, update LM state, recombine tokens
*/
void Decoder::HandleWordend (LexNode *ln)
{
   LexNodeInst *inst;
   WordendHyp *weHyp, *prev;
   TokenSet *ts;
   RelToken *tok, *tokJ;
   int i, j, newN;
   LMState dest;
   LMTokScore lmScore;
   PronId pronid;
   RelTokScore newDelta, bestDelta, deltaLimit;
#ifdef MODALIGN
   ModendHyp *modpath;
#endif

   assert (ln->type == LN_WORDEND);

   inst = ln->inst;
   ts = inst->ts;

   /* main beam pruning for reltoks */
   /* main and relative beam pruning */
   deltaLimit = std::max(_decInst->beamLimit - ts->score, _decInst->relBeamWidth);

   /* for each token i in set, take transition in LM
      recombine tokens in same LMState 
      newN is (current) number of new tokens (newN <= ts->n) */

   pronid = (PronId) ln->data.pron;

   newN = 0;
   bestDelta = LZERO;
   for (i = 0; i < ts->n; ++i) {
      tok = &ts->relTok[i];

      if (tok->delta < deltaLimit)
         continue;      /* prune */

      lmScore = LMCacheTransProb (_decInst->lm, tok->lmState, pronid, &dest);

      /* word insertion penalty */
      lmScore += _decInst->insPen;
      /* remember prev path now, as we might overwrite it below */
      prev = tok->path;
#ifdef MODALIGN
      modpath = tok->modpath;
#endif

      /* subtract lookahead which has already been applied */
      newDelta = tok->delta + (lmScore - tok->lmscore);

      if (newDelta < deltaLimit)
         continue;              /* prune */

      bestDelta = std::max(bestDelta, newDelta);

      /* insert in list */
      for (j = 0; j < newN; ++j) {      /* is there already a token in state dest? */
         tokJ = &ts->relTok[j];
         if (tokJ->lmState == dest) {
            if (!_decInst->latgen) {
               if (newDelta > tokJ->delta) {        /* replace tokJ */
                  tokJ->delta = newDelta;
                  tokJ->lmscore = 0.0;     /* reset lookahead */
#ifdef MODALIGN
                  tokJ->modpath = modpath;
#endif
                  /* update path; 
                     weHyp exists, pron is the same anyway, update rest */
                  tokJ->path->prev = prev;
                  tokJ->path->score = ts->score + newDelta;
                  tokJ->path->lm = lmScore;
#ifdef MODALIGN
                  tokJ->path->modpath = modpath;
#endif
               }
               /* else just toss token */
            }
            else {      /* latgen */
               AltWordendHyp *alt;

               alt = (AltWordendHyp *) New (&_decInst->altweHypHeap, sizeof (AltWordendHyp));
               
               if (newDelta > tokJ->delta) {
                  /* move tokJ->path to alt */
                  alt->prev = tokJ->path->prev;
                  alt->score = tokJ->path->score;
                  alt->lm = tokJ->path->lm;
#ifdef MODALIGN
                  alt->modpath = tokJ->path->modpath;
#endif

                  /* replace tokJ */
                  tokJ->delta = newDelta;
                  tokJ->lmscore = 0.0;     /* reset lookahead */
#ifdef MODALIGN
                  tokJ->modpath = modpath;

                  tokJ->path->modpath = modpath;
#endif

                  /* store new tok info in path 
                     weHyp exists, pron is the same anyway, update rest */
                  tokJ->path->prev = prev;
                  tokJ->path->score = ts->score + newDelta;
                  tokJ->path->lm = lmScore;
               }
               else {
                  /* store new tok info in alt */
                  alt->prev = prev;
                  alt->score = ts->score + newDelta;
                  alt->lm = lmScore;
#ifdef MODALIGN
                  alt->modpath = modpath;
#endif
               }

               /* attach alt to tokJ's weHyp */
               alt->next = tokJ->path->alt;
               tokJ->path->alt = alt;
            }               
            break;      /* leave j loop */
         }
      }

      if (j == newN) {          /* no token in state dest yet */
         int k;

         /* find spot to insert LMState dest */
         for (j = 0; j < newN; ++j)
            if (ts->relTok[j].lmState > dest)
               break;

         /* move any following reltokens up one slot */
         for (k = newN ; k > j; --k)
            ts->relTok[k] = ts->relTok[k-1];
         
         tokJ = &ts->relTok[j];
         ++newN;

         /* new wordendHyp */
         weHyp = (WordendHyp *) New (&_decInst->weHypHeap, sizeof (WordendHyp));
      
         weHyp->prev = prev;
         weHyp->pron = ln->data.pron;
         weHyp->score = ts->score + newDelta;
         weHyp->lm = lmScore;
         weHyp->frame = _decInst->frame;
         weHyp->alt = NULL;
         weHyp->user = 0;
#ifdef MODALIGN
         weHyp->modpath = modpath;

         tokJ->modpath = modpath;
#endif

         tokJ->path = weHyp;
         /* only really necessary if (i!=j) i.e. (tok!=tokJ) */
         tokJ->delta = newDelta;
         tokJ->lmState = dest;
         tokJ->we_tag = (void *) ln;
         tokJ->lmscore = 0.0;   /* reset lookahead */
      }
   } /* for token i */

   ts->n = newN;

   if (newN > 0) {
      AltWordendHyp *alt;
      /* renormalise  to new best score */
      for (i = 0; i < ts->n; ++i) {
         tok = &ts->relTok[i];
         tok->delta -= bestDelta;

         /* convert alt wordendHyp scores to deltas relativ to main weHyp */
         for (alt = tok->path->alt; alt; alt = alt->next)
            alt->score = alt->score - tok->path->score;
      }
      ts->score += bestDelta;

      /* ####  TokSet id */
      ts->id = ++_decInst->tokSetIdCount;
   }
   else {
      ts->id = 0;
      ts->score = LZERO;
   }

   inst->best = ts->score;
}


/* UpdateWordEndHyp

     update wordend hyps of all tokens with current time and score
 */
void Decoder::UpdateWordEndHyp ( LexNodeInst *inst)
{
   int i;
   TokenSet *ts;
   RelToken *tok;
   WordendHyp *weHyp, *oldweHyp;
   
   ts = inst->ts;
           
   for (i = 0; i < ts->n; ++i) {
      tok = &ts->relTok[i];

      oldweHyp = tok->path;

      /* don't copy weHyp, if it is up-to-date (i.e. for <s>) */
      if (oldweHyp->frame != _decInst->frame || oldweHyp->pron != _decInst->net->startPron) {
         weHyp = (WordendHyp *) New (&_decInst->weHypHeap, sizeof (WordendHyp));
         *weHyp = *oldweHyp;
         weHyp->score = ts->score + tok->delta;
         weHyp->frame = _decInst->frame;
#ifdef MODALIGN
         weHyp->modpath = tok->modpath;
#endif
         tok->path = weHyp;

         /* altweHyps don't need to be changed  here as they are relative to 
            the main weHyp's score*/
      }
#ifdef MODALIGN
      tok->modpath = NULL;
#endif
   }
}

void Decoder::AddPronProbs (TokenSet *ts, int var)
{
   /* #### this is pathetically slow!!
      #### fix by storing 3 lists parallel to pronlist with 
      #### prescaled pronprobs
   */
   RelTokScore bestDelta = LZERO;

   for (int i = 0; i < ts->n; ++i) {
      RelToken* tok = &(ts->relTok[i]);
      WordendHyp *path = tok->path;
      Pron pron = _decInst->net->pronlist[path->pron];

      if (var == 1)
         pron = pron->next;
      else if (var == 2)
         pron = pron->next->next;
      
      tok->delta += _decInst->pronScale * pron->prob;
      if (tok->delta > bestDelta)
         bestDelta = tok->delta;

      /* need to make copy of path before modifying it */
      if (path->user != var) {
         WordendHyp *weHyp;

         weHyp = (WordendHyp *) New (&_decInst->weHypHeap, sizeof (WordendHyp));
         *weHyp = *path;
         weHyp->user = var;
         tok->path = weHyp;
      }
   }

   /* renormalise token set */
   for (int i = 0; i < ts->n; ++i)
      ts->relTok[i].delta -= bestDelta;
   ts->score += bestDelta;
}


void Decoder::HandleSpSkipLayer (LexNodeInst *inst)
{
   LexNode *ln;

   ln = inst->node;
   if (ln->nfoll == 1) {    /* sp is unique foll, for sil case there are two! */
      LexNode *lnSA;
      
      assert (ln->foll[0]->data.hmm == _decInst->net->hmmSP);
      assert (ln->foll[0]->nfoll == 1);

      /* bypass sp model for - variant */
      /*  propagate tokens to follower of sp (ln->foll[0]->foll[0]) */
      /*    adding - variant pronprob (pron->prob) */
      lnSA = ln->foll[0]->foll[0];
                  
      /* node should be either inactive or empty */
      assert (!lnSA->inst || lnSA->inst->ts[0].n == 0);
      
      PropIntoNode (&inst->ts[0], ln->foll[0]->foll[0], FALSE);
      
      /* add pronprobs and keep record of variant in path->user */
      /*   user = 0: - variant, 1: sp, 2: sil */
      AddPronProbs (&lnSA->inst->ts[0], 0);
      
      /* now add sp variant pronprob to token set and propagate as normal */
      AddPronProbs (&inst->ts[0], 1);
      PropagateExternal (inst, FALSE, FALSE);
   }
   else {   /* sil variant */
      int sentEnd = 0;
      TokenSet *tempTS;

      tempTS = _decInst->tempTS[1];
      tempTS->score = 0.0;
      tempTS->n = 0;
      tempTS->id = 0;

      assert (ln->nfoll == 2);  /* inter word sil and path to sent end */
      
      if (ln->foll[1]->type == LN_CON)  /* lnTime node in SA, see comment 
                                           in HLVNet:CreateStartEnd() */
         sentEnd = 1;
      
      /*   path to SENT_END */
      /* propagate to ln->foll[sentEnd] and add - var pronpob */
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 0);
      if (tempTS->score >= _decInst->beamLimit)
         PropIntoNode (tempTS, ln->foll[sentEnd], FALSE);
      /* propagate to SENT_END  sp */
      tempTS->score = 0.0; tempTS->n = 0; tempTS->id = 0;
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 1);
      if (tempTS->score >= _decInst->beamLimit)
         PropIntoNode (tempTS, _decInst->net->lnSEsp, FALSE);
      /* propagate to SENT_END  sil */
      tempTS->score = 0.0; tempTS->n = 0; tempTS->id = 0;
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 2);
      if (tempTS->score >= _decInst->beamLimit)
         PropIntoNode (tempTS, _decInst->net->lnSEsil, FALSE);
      
      
      /*   normal word loop */
      /* add sil variant pronprob to token set and propagate */
      AddPronProbs (&inst->ts[0], 2);
      if (inst->ts[0].score < _decInst->beamLimit) { /* prune to keep MTS happy */
         inst->ts[0].n = 0;
         inst->ts[0].id = 0;
      }
      else {
         PropIntoNode (&inst->ts[0], ln->foll[1 - sentEnd], FALSE);
      }
   }
}

void Decoder::ZSLayerBeamPruning(LexNodeInst * head, TokScore& beamLimit) {
  TokScore bestScore = LZERO;
  for (; head; head = head->next )
    bestScore = std::max(bestScore, head->best);

  beamLimit = std::max(bestScore - _decInst->zsBeamWidth, _decInst->beamLimit);
}


inline void eraseLinkedListNode(LexNodeInst* &head, LexNodeInst* prev, LexNodeInst* inst) {
  if (prev)
    prev->next = inst->next;
  else                 // first inst in layer
    head = inst->next;
}

void Decoder::PropagateInternal() {
   /* internal token propagation:
      order doesn't really matter, but we use the same as for external propagation */
   for (int l = 0; l < _decInst->nLayers; ++l) {
      for (LexNodeInst* inst = _decInst->instsLayer[l]; inst; inst = inst->next) {
	 switch (inst->node->type) {
	   case LN_MODEL:   /* Model node */
	     PropagateInternal (inst);
	     break;
	   case LN_CON:	    /* Context or Wordend node */
	   case LN_WORDEND:
	     /* clear tokenset in preparation for external propagation*/
	     inst->ts[0].n = 0;
	     inst->best = LZERO;
	     break;
	 }
      }
   }
}

void Decoder::SetObservation(Observation **obsBlock, int nObs) {
  _decInst->obs = obsBlock[0];
  _decInst->nObs = nObs;
  for (int i = 0; i < nObs; ++i)
    _decInst->obsBlock[i] = obsBlock[i];
  _decInst->bestScore = LZERO;
  _decInst->bestInst = NULL;
  ++_decInst->frame;
}

void Decoder::WordEndBeamPruning(LexNodeInst* head, TokScore &beamLimit) {
  TokScore bestWEscore = LZERO;
  LexNodeInst *next;

  LexNodeInst* prev = NULL;
  for (auto inst = head; inst; inst = next) {
    next = inst->next;     /* store now, we might free inst below! */

    if (inst->best < beamLimit) {  /* global main beam */
      eraseLinkedListNode(head, prev, inst);
      DeactivateNode (inst->node);
    }
    else {
      HandleWordend (inst->node);
      bestWEscore = std::max(bestWEscore, inst->best);
      prev = inst;
    }
  }

  beamLimit = std::max(bestWEscore - _decInst->weBeamWidth, _decInst->beamLimit); /* global beam is tighter */
}

void Decoder::RelaxBeamLimit() {

   const size_t MMP_NBINS = 128;
   if (_decInst->maxModel > 0) {
      int hist[MMP_NBINS];

      for (int i = 0; i < MMP_NBINS; ++i)
         hist[i] = 0;
      
      int nModel = 0;
      LogFloat binWidth = _decInst->curBeamWidth / MMP_NBINS;
      /* fill histogram */
      for (int l = 0; l < _decInst->nLayers; ++l) {
         for (auto inst = _decInst->instsLayer[l]; inst; inst = inst->next) {
            if (inst->best <= LSMALL)
	      continue;

	    int bin = (_decInst->bestScore - inst->best) / binWidth;
	    if (bin < MMP_NBINS) {
	      ++hist[bin];
	      ++nModel;
	    }
         }
      }
      
      if (nModel <= _decInst->maxModel) {
         /* slowly increase beamWidth again */
         _decInst->curBeamWidth *= dynBeamInc;
         if (_decInst->curBeamWidth > _decInst->beamWidth)
            _decInst->curBeamWidth = _decInst->beamWidth;
      }
   }

   _decInst->beamLimit = _decInst->bestScore - _decInst->curBeamWidth;
}

/* ProcessFrame

     Takes the observation vector and propatagets all tokens and
     performs pruning as necessary.
*/
void Decoder::ProcessFrame (Observation **obsBlock, int nObs, AdaptXForm *xform) {
   TokScore beamLimit;
   
   inXForm = xform; /* sepcifies the transform to use */
   
   /* reset obs */
   SetObservation(obsBlock, nObs);

   GarbageCollectPaths ();

   PropagateInternal();

   /* now for all LN_MODEL nodes inst->best is set, this is used to determine 
      the lower beam limit */

   _decInst->beamLimit = _decInst->bestScore - _decInst->curBeamWidth;

   /* beam pruning & external propagation */
   for (int l = 0; l < _decInst->nLayers; ++l) {
      LexNodeInst* &head = _decInst->instsLayer[l];

      /* update word end time and score in tok->path when passing
         through the appropriate layer */
      if (l == _decInst->net->wordEndLayerId) {
         for (auto inst = head; inst; inst = inst->next)
            UpdateWordEndHyp (inst);
      }

      beamLimit = _decInst->beamLimit;
      if ((_decInst->weBeamWidth < _decInst->beamWidth) && (l == LAYER_WE))
	WordEndBeamPruning(head, beamLimit);
      else if ((_decInst->zsBeamWidth < _decInst->beamWidth) && (l == LAYER_ZS || l == LAYER_SA))
	ZSLayerBeamPruning(head, beamLimit);

      /* Due to the layer-by-layer structure inst->best values for
         LN_CON and LN_WORDEND nodes are set before they are examined
         for pruning purposes, although they were not available when the 
         beamlimit  was set.
      */
      LexNodeInst* prev = NULL, *next = NULL;
      for (auto inst = head; inst; inst = next) {
	 /* store now, we might free inst below! */
         next = inst->next;     
         
         if (inst->node->type != LN_WORDEND && inst->node->lmlaIdx != 0 && inst->ts->n > 0) {
	    /* don't bother if inst will be pruned anyway */
            if (inst->ts[0].score >= beamLimit)       
               UpdateLMlookahead (inst->node);

	    /* UpLMLA might have killed the entire TS, esp. in latlm */
            if (inst->ts->n > 0)
	      inst->best = std::max(inst->best, inst->ts[0].score);
         }

         if (inst->best < beamLimit) { /* take inst out of instsLayer list and deactivate it */
	   eraseLinkedListNode(head, prev, inst);
	   DeactivateNode (inst->node);
	   continue;
         }

#ifdef COLLECT_STATS
	 ++_decInst->stats.nActive;
#endif
	 /* special code for pronprob handling before sil layer */
	 if (_decInst->net->silDict && (l == _decInst->net->spSkipLayer) && inst->ts[0].n > 0) {
	   HandleSpSkipLayer (inst);
	   continue;
	 }
	 
	 /* normal case: non silDict or non spSkipLayer */
	 /* call HandleWordend, if we don't do we-pruning or we are in in LAYER_SIL and LAYER_AB (where we need to handle the wordend nodes for SENT_START and SENT_END). */
	 /* ### fix this */
	 /* experiment for richer lattices. keep sp and sil variants distinct by marking sil in LSBit of tok->we_tag */ /* #### we need the equivalent for pronprob sildicts! */
	 if (l == LAYER_SIL && inst->node->type == LN_MODEL) {
	   auto &H = inst->node->data.hmm;
	   if (H != _decInst->net->hmmSP) {
	     int N = H->numStates;
	     TokenSet *ts = &inst->ts[N-1];
	     for (int i = 0; i < ts->n; ++i)
	       ts->relTok[i].set_we_tag();
	   }
	 }

	 bool handleWE = !(_decInst->weBeamWidth < _decInst->beamWidth) || (l == LAYER_SIL) || (l == LAYER_AB);
	 PropagateExternal (inst, handleWE, l == LAYER_BY);

	 prev = inst;
      } /* for inst */
   }    /* for layer */
  
   this->RelaxBeamLimit();

#ifdef COLLECT_STATS
   ++_decInst->stats.nFrames;
#endif

   PI_LR = PI_GEN = 0;
   _decInst->outPCache->cacheHit = _decInst->outPCache->cacheMiss = 0;
   _decInst->lmCache->transHit = _decInst->lmCache->transMiss = 0;
   _decInst->lmCache->laHit = _decInst->lmCache->laMiss = 0;

   if (_decInst->nPhone > 0)
      CalcPhonePost ();
}


/*  CC-mode style info for emacs
 Local Variables:
 c-file-style: "htk"
 End:
*/
