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

/* perform Bucket sort/Histogram pruning to reduce to _dec->nTok tokens */
void Decoder::TokSetHistogramPruning(TokenSet *dest, const RelTokScore& deltaLimit,
    int nWinTok, RelToken* &winTok, TokScore &winScore) {

  const size_t NBINS = 64;


  dest->id = ++_dec->tokSetIdCount;    /* #### new id always necessary? */
  LogFloat binWidth = deltaLimit * 1.001 / NBINS;   /* handle delta == deltaLimit case */

  int n[NBINS] = {0};

  for (int i = 0; i < nWinTok; ++i)
    ++n[(int) (winTok[i].delta / binWidth)];

  int nTok = 0;
  int i = -1, j;
  while (nTok < _dec->nTok)
    nTok += n[++i];

  if (nTok == _dec->nTok) {
    int binLimit = i;

    for (i = 0, j = 0; i < nWinTok; ++i) { 
      if ((int) (winTok[i].delta / binWidth) <= binLimit) {
	dest->relTok[j] = winTok[i];
	++j;
      }
    }
  }
  else {

    /* do not include last bin */
    LogFloat limit = binWidth * i;
    nTok -= n[i]; 

    /* need to relax limit so that we get an extra (_dec->nTok - nTok) tokens */
    /* #### very simplistic implementation -- imporve? */

    LogFloat bestDelta = limit;
    int nBetter;
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
    } while (nBetter < _dec->nTok);

    /* multiple tokens with delta == limit ==> delete some */
    if (nBetter > _dec->nTok) {  
      for (i = 0; nBetter > _dec->nTok; ++i)
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

  dest->n = _dec->nTok;
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

void Decoder::CopyWinningToken(TokenSet *src, TokenSet *dest, int nWinTok, RelToken *winTok,
    const TokScore &winScore, const int *nWin) {

     /* just copy */
     for (int i = 0; i < nWinTok; ++i)
       dest->relTok[i] = winTok[i];

     dest->n = nWinTok;
     dest->score = winScore;

     if (nWinTok == nWin[0])
       dest->id = src->id;          /* copy src->id */
     else if (nWinTok == nWin[1])
       dest->id = dest->id;         /* copy dest->id */
     else
       dest->id = ++_dec->tokSetIdCount;    /* new id */
}

/* MergeTokSet

     Merge TokenSet src into dest after adding score to all src scores
     by recombining tokens in the same LM state and keeping only the
     _dec->nTok best tokens in different LM states.

*/
void Decoder::MergeTokSet (TokenSet *src, TokenSet *dest, LogFloat score, bool prune) {
   /* empty dest tokenset -> just copy src */
   if (dest->n == 0) {   
      *dest = *src;
      dest->score += score;
      return;

   } else if (prune && src->score + score < _dec->beamLimit)
     return;

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
   deltaLimit = _dec->nTok * _dec->relBeamWidth;  

   /* set pruning deltaLimit */
   if (prune)
     deltaLimit = std::max(_dec->beamLimit - winScore, _dec->relBeamWidth);

   /* find winning tokens */
   /* location where winnning toks came from: 0 == src, 1 == dest */
   RelToken *winTok = _dec->winTok;
   int nWinTok = 0;
   int nWin[2] = {0, 0};
   FindWinningTokens(src, dest, score, srcCorr, destCorr, deltaLimit, winTok, nWinTok, nWin);

   if (nWinTok <= _dec->nTok)
     CopyWinningToken(src, dest, nWinTok, winTok, winScore, nWin);
   else 
     TokSetHistogramPruning(dest, deltaLimit, nWinTok, winTok, winScore);
}

void Decoder::__collect_stats__ (TokenSet *instTS, int N) {
  for (int i = 0; i < N; ++i) {
    if (instTS[i].n > 0) {
      ++_dec->stats.nTokSet;
      _dec->stats.sumTokPerTS += instTS[i].n;
    }
  }
}

void Decoder::OptLeftToRightPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm) {

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
      LogFloat outP = cOutP (_dec->obs, hmm, i);
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

  if (bestScore > _dec->bestScore) {
    _dec->bestScore = bestScore;
    _dec->bestInst = inst;
  }

  __collect_stats__(instTS, N);
}

void Decoder::GeneralPropagateInternal(LexNodeInst *inst, TokenSet* instTS, int N, SMatrix &trP, HLink &hmm) {

  /* temp storage for N tokensets  */
  TokenSet *tempTS = _dec->tempTS[N];

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
      LogFloat outP = cOutP (_dec->obs, hmm, j);
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
  if (bestScore > _dec->bestScore) {
    _dec->bestScore = bestScore;
    _dec->bestInst = inst;
  }

  /* # this only collects stats for the model nodes */
  __collect_stats__(instTS, N);
}

/* \brief InternalPropagation
 *    Internal token propagation
 *    ln->type must be LN_MODEL
 * */
void Decoder::InternalPropagation (LexNodeInst *inst) {

   LexNode *ln = inst->node;
   HLink hmm = ln->data.hmm;
   int N = hmm->numStates;
   SMatrix trP = hmm->transP;
   TokenSet* instTS = inst->ts;

   /* Main beam pruning: prune tokensets before propagation
      the beamLimit is the one found during the last frame */
   auto ts = &instTS[0];
   for (int i = 1; i < N; ++i, ++ts) {
     if (ts->score < _dec->beamLimit)
       ts->n = ts->id = 0;
   }

   /* Optimised version for L-R models OR 
      General (not left-to-right) Propagate Internal   */
   if (hmm->tIdx < 0)
      OptLeftToRightPropagateInternal(inst, instTS, N, trP, hmm);
   else
     GeneralPropagateInternal(inst, instTS, N, trP, hmm);
}

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
      m = (ModendHyp*) New (&_dec->modendHypHeap, sizeof (ModendHyp));
      m->frame = _dec->frame;
      m->ln = ln;
      m->prev = tok->modpath;

      tok->modpath = m;
   }
}

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
          (entryTS->score > _dec->beamLimit)) {
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
      if (exitTS->score < _dec->beamLimit) {
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
   if (exitTS->n > 0 && exitTS->score > _dec->beamLimit) {
      if (_dec->modAlign)
         UpdateModPaths (exitTS, ln);
      /* loop over following nodes */
      for (i = 0; i < ln -> nfoll; ++i) {
         follLN = ln->foll[i];
         PropIntoNode (exitTS, follLN, wintTree);
      } /* for i following nodes */
   }
}


/* HandleWordend
     update traceback, add LM, update LM state, recombine tokens
     ln->type must be LN_WORDEND
*/
void Decoder::HandleWordend (LexNode *ln)
{
   WordendHyp *prev;
   LMState dest;
   PronId pronid;
   RelTokScore newDelta, bestDelta, deltaLimit;
   ModendHyp *modpath;

   LexNodeInst *inst = ln->inst;
   TokenSet *ts = inst->ts;

   /* main beam pruning for reltoks */
   /* main and relative beam pruning */
   deltaLimit = std::max(_dec->beamLimit - ts->score, _dec->relBeamWidth);

   /* for each token i in set, take transition in LM
      recombine tokens in same LMState 
      newN is (current) number of new tokens (newN <= ts->n) */

   pronid = (PronId) ln->data.pron;

   int newN = 0;
   bestDelta = LZERO;
   for (int i = 0; i < ts->n; ++i) {
      auto tok = &ts->relTok[i];

      // Pruning
      if (tok->delta < deltaLimit) continue;      

      auto lmScore = LMCacheTransProb (_dec->lm, tok->lmState, pronid, &dest);

      // word insertion penalty
      lmScore += _dec->insPen;

      /* remember prev path now, as we might overwrite it below */
      prev = tok->path;
      modpath = tok->modpath;

      /* subtract lookahead which has already been applied */
      newDelta = tok->delta + (lmScore - tok->lmscore);

      // Prune again
      if (newDelta < deltaLimit) continue;

      bestDelta = std::max(bestDelta, newDelta);

      int j;
      /* insert in list */
      for (j = 0; j < newN; ++j) {      /* is there already a token in state dest? */
         auto tokJ = &ts->relTok[j];

         if (tokJ->lmState != dest)
	   continue;

	 if (!_dec->latgen) {
	   // replace tokJ
	   if (newDelta > tokJ->delta) {
	     tokJ->delta = newDelta;
	     tokJ->lmscore = 0.0;	  // reset lookahead
	     tokJ->modpath = modpath;

	     // Update path -- weHyp exists, pron is the same anyway, update rest
	     tokJ->path->prev = prev;
	     tokJ->path->score = ts->score + newDelta;
	     tokJ->path->lm = lmScore;
	     tokJ->path->modpath = modpath;
	   }
	   /* else just toss token */
	 }
	 else {      /* latgen */
	   auto alt = (AltWordendHyp *) New (&_dec->altweHypHeap, sizeof (AltWordendHyp));

	   if (newDelta > tokJ->delta) {
	     /* move tokJ->path to alt */
	     alt->prev = tokJ->path->prev;
	     alt->score = tokJ->path->score;
	     alt->lm = tokJ->path->lm;
	     alt->modpath = tokJ->path->modpath;

	     /* replace tokJ */
	     tokJ->delta = newDelta;
	     tokJ->lmscore = 0.0;     /* reset lookahead */
	     tokJ->modpath = modpath;

	     tokJ->path->modpath = modpath;

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
	     alt->modpath = modpath;
	   }

	   /* attach alt to tokJ's weHyp */
	   alt->next = tokJ->path->alt;
	   tokJ->path->alt = alt;
	 }               
	 break;      /* leave j loop */
      }

      if (j == newN) {          /* no token in state dest yet */
         /* find spot to insert LMState dest */
         for (j = 0; j < newN; ++j)
            if (ts->relTok[j].lmState > dest)
               break;

         /* move any following reltokens up one slot */
         for (int k = newN ; k > j; --k)
            ts->relTok[k] = ts->relTok[k-1];
         
         auto tokJ = &ts->relTok[j];
         ++newN;

         /* new wordendHyp */
         auto weHyp = (WordendHyp *) New (&_dec->weHypHeap, sizeof (WordendHyp));
      
         weHyp->prev = prev;
         weHyp->pron = ln->data.pron;
         weHyp->score = ts->score + newDelta;
         weHyp->lm = lmScore;
         weHyp->frame = _dec->frame;
         weHyp->alt = NULL;
         weHyp->user = 0;
         weHyp->modpath = modpath;

         tokJ->modpath = modpath;
         tokJ->path = weHyp;
         tokJ->delta = newDelta;
         tokJ->lmState = dest;
         tokJ->we_tag = (void *) ln;
         tokJ->lmscore = 0.0;
      }
   } /* for token i */

   ts->n = newN;

   if (newN > 0) {
      /* renormalise  to new best score */
      for (int i = 0; i < ts->n; ++i) {
         auto tok = &ts->relTok[i];
         tok->delta -= bestDelta;

         /* convert alt wordendHyp scores to deltas relativ to main weHyp */
         for (auto alt = tok->path->alt; alt; alt = alt->next)
            alt->score = alt->score - tok->path->score;
      }
      ts->score += bestDelta;

      ts->id = ++_dec->tokSetIdCount;
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
void Decoder::UpdateWordEndHyp ( LexNodeInst *inst) {

   auto *ts = inst->ts;
           
   for (int i = 0; i < ts->n; ++i) {
      auto tok = &ts->relTok[i];
      auto oldweHyp = tok->path;

      /* don't copy weHyp, if it is up-to-date (i.e. for <s>) */
      if (oldweHyp->frame != _dec->frame || oldweHyp->pron != _dec->net->startPron) {
         auto weHyp = (WordendHyp *) New (&_dec->weHypHeap, sizeof (WordendHyp));
         *weHyp = *oldweHyp;
         weHyp->score = ts->score + tok->delta;
         weHyp->frame = _dec->frame;
         weHyp->modpath = tok->modpath;
         tok->path = weHyp;
      }
      tok->modpath = NULL;
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
      Pron pron = _dec->net->pronlist[path->pron];

      if (var == 1)
         pron = pron->next;
      else if (var == 2)
         pron = pron->next->next;
      
      tok->delta += _dec->pronScale * pron->prob;
      if (tok->delta > bestDelta)
         bestDelta = tok->delta;

      /* need to make copy of path before modifying it */
      if (path->user != var) {
         WordendHyp *weHyp;

         weHyp = (WordendHyp *) New (&_dec->weHypHeap, sizeof (WordendHyp));
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

/* \brief HandleSpSkipLayer
 * Handle short pause (sp) and skip layer
 *
 * */
void Decoder::HandleSpSkipLayer (LexNodeInst *inst)
{
   LexNode *ln = inst->node;

   if (ln->nfoll == 1) {
      LexNode *lnSA = ln->foll[0]->foll[0];
                  
      PropIntoNode (&inst->ts[0], ln->foll[0]->foll[0], FALSE);
      
      AddPronProbs (&lnSA->inst->ts[0], 0);
      AddPronProbs (&inst->ts[0], 1);

      PropagateExternal (inst, FALSE, FALSE);
   }
   else {   /* sil variant */
      int sentEnd = 0;

      TokenSet *tempTS = _dec->tempTS[1];
      tempTS->score = tempTS->n = tempTS->id = 0;

      /* lnTime node in SA, see comment in HLVNet:CreateStartEnd() */
      if (ln->foll[1]->type == LN_CON)  
         sentEnd = 1;
      
      /*   path to SENT_END */
      /* propagate to ln->foll[sentEnd] and add - var pronpob */
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 0);
      if (tempTS->score >= _dec->beamLimit)
         PropIntoNode (tempTS, ln->foll[sentEnd], FALSE);

      /* propagate to SENT_END  sp */
      tempTS->score = 0.0; tempTS->n = 0; tempTS->id = 0;
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 1);
      if (tempTS->score >= _dec->beamLimit)
         PropIntoNode (tempTS, _dec->net->lnSEsp, FALSE);

      /* propagate to SENT_END  sil */
      tempTS->score = 0.0; tempTS->n = 0; tempTS->id = 0;
      MergeTokSet (&inst->ts[0], tempTS, 0.0, FALSE);
      AddPronProbs (tempTS, 2);
      if (tempTS->score >= _dec->beamLimit)
         PropIntoNode (tempTS, _dec->net->lnSEsil, FALSE);
      
      // Normal word loop -- add sil variant pronprob to token set and propagate
      AddPronProbs (&inst->ts[0], 2);
      if (inst->ts[0].score < _dec->beamLimit)
         inst->ts[0].n = inst->ts[0].id = 0;
      else
         PropIntoNode (&inst->ts[0], ln->foll[1 - sentEnd], FALSE);
   }
}

void Decoder::ZSLayerBeamPruning(LexNodeInst * head, TokScore& beamLimit) {
  TokScore bestScore = LZERO;
  for (; head; head = head->next )
    bestScore = std::max(bestScore, head->best);

  beamLimit = std::max(bestScore - _dec->zsBeamWidth, _dec->beamLimit);
}

inline void RemoveLexNode(LexNodeInst* &head, LexNodeInst* prev, LexNodeInst* inst) {
  if (prev)
    prev->next = inst->next;
  else
    head = inst->next;
}

/* \brief InternalPropagation
 *
 * */

void Decoder::InternalPropagation() {
  // static int counter = 0; if (++counter % 20 == 1) printf("\33[33mLAYER_Z   LAYER_ZS  LAYER_SIL LAYER_SA  LAYER_A   LAYER_AB  LAYER_BY  LAYER_WE  LAYER_YZ\33[0m\n");

   for (int l = 0; l < _dec->nLayers; ++l) {
     int nActive = 0;
      for (LexNodeInst* inst = _dec->instsLayer[l]; inst; inst = inst->next) {
	++nActive;
	 switch (inst->node->type) {
	   case LN_MODEL:   /* Model node */
	     InternalPropagation (inst);
	     break;
	   case LN_CON:	    /* Context or Wordend node */
	   case LN_WORDEND:
	     /* clear tokenset in preparation for external propagation*/
	     inst->ts[0].n = 0;
	     inst->best = LZERO;
	     break;
	 }
      }
      // printf("%5d     ", nActive);
   }
   // printf("\33[0m \n");
}

void Decoder::SetObservation(Observation **obsBlock, int nObs) {
  _dec->obs = obsBlock[0];
  _dec->nObs = nObs;
  for (int i = 0; i < nObs; ++i)
    _dec->obsBlock[i] = obsBlock[i];
  _dec->bestScore = LZERO;
  _dec->bestInst = NULL;
  ++_dec->frame;
}

void Decoder::WordEndBeamPruning(LexNodeInst* head, TokScore &beamLimit) {
  TokScore bestWEscore = LZERO;
  LexNodeInst *next;

  LexNodeInst* prev = NULL;
  for (auto inst = head; inst; inst = next) {
    next = inst->next;     /* store now, we might free inst below! */

    if (inst->best < beamLimit) {  /* global main beam */
      RemoveLexNode(head, prev, inst);
      DeactivateNode (inst->node);
    }
    else {
      HandleWordend (inst->node);
      bestWEscore = std::max(bestWEscore, inst->best);
      prev = inst;
    }
  }

  beamLimit = std::max(bestWEscore - _dec->weBeamWidth, _dec->beamLimit); /* global beam is tighter */
}

void Decoder::RelaxBeamLimit() {

   const size_t MMP_NBINS = 128;
   if (_dec->maxModel > 0) {
      int hist[MMP_NBINS] = {0};

      int nModel = 0;
      LogFloat binWidth = _dec->curBeamWidth / MMP_NBINS;
      /* fill histogram */
      for (int l = 0; l < _dec->nLayers; ++l) {
         for (auto inst = _dec->instsLayer[l]; inst; inst = inst->next) {
            if (inst->best <= LSMALL)
	      continue;

	    int bin = (_dec->bestScore - inst->best) / binWidth;
	    if (bin < MMP_NBINS) {
	      ++hist[bin];
	      ++nModel;
	    }
         }
      }
      
      if (nModel <= _dec->maxModel) {
         /* slowly increase beamWidth again */
         _dec->curBeamWidth *= dynBeamInc;
         if (_dec->curBeamWidth > _dec->beamWidth)
            _dec->curBeamWidth = _dec->beamWidth;
      }
   }

   _dec->beamLimit = _dec->bestScore - _dec->curBeamWidth;
}

void Decoder::SetWordEndTag(LexNodeInst* inst) {

  auto &H = inst->node->data.hmm;
  if (H != _dec->net->hmmSP) {
    int N = H->numStates;
    TokenSet *ts = &inst->ts[N-1];
    for (int i = 0; i < ts->n; ++i)
      ts->relTok[i].set_we_tag();
  }

}

/* ProcessFrame

     Takes the observation vector and propatagets all tokens and
     performs pruning as necessary.
*/
void Decoder::ProcessFrame (Observation **obsBlock, int nObs, AdaptXForm *xform) {
   
   inXForm = xform; /* sepcifies the transform to use */
   
   /* reset obs */
   SetObservation(obsBlock, nObs);

   GarbageCollectPaths ();

   InternalPropagation();

   /* now for all LN_MODEL nodes inst->best is set, this is used to determine 
      the lower beam limit */

   _dec->beamLimit = _dec->bestScore - _dec->curBeamWidth;

   /* beam pruning & external propagation */
   for (int l = 0; l < _dec->nLayers; ++l) {
      LexNodeInst* &head = _dec->instsLayer[l];

      /* update word end time and score in tok->path when passing
         through the appropriate layer */
      if (l == _dec->net->wordEndLayerId) {
         for (auto inst = head; inst; inst = inst->next)
            UpdateWordEndHyp (inst);
      }

      TokScore beamLimit = _dec->beamLimit;
      if ((_dec->weBeamWidth < _dec->beamWidth) && (l == LAYER_WE))
	WordEndBeamPruning(head, beamLimit);
      else if ((_dec->zsBeamWidth < _dec->beamWidth) && (l == LAYER_ZS || l == LAYER_SA))
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
	   RemoveLexNode(head, prev, inst);
	   DeactivateNode (inst->node);
	   continue;
         }

	 ++_dec->stats.nActive;
	 /* special code for pronprob handling before sil layer */
	 if (_dec->net->silDict && (l == _dec->net->spSkipLayer) && inst->ts[0].n > 0) {
	   HandleSpSkipLayer (inst);
	   continue;
	 }
	 
	 /* normal case: non silDict or non spSkipLayer */
	 /* call HandleWordend, if we don't do we-pruning or we are in in LAYER_SIL and LAYER_AB (where we need to handle the wordend nodes for SENT_START and SENT_END). */
	 /* ### fix this */
	 /* experiment for richer lattices. keep sp and sil variants distinct by marking sil in LSBit of tok->we_tag */ /* #### we need the equivalent for pronprob sildicts! */
	 if (l == LAYER_SIL && inst->node->type == LN_MODEL)
	   this->SetWordEndTag(inst);

	 bool handleWE = !(_dec->weBeamWidth < _dec->beamWidth) || (l == LAYER_SIL) || (l == LAYER_AB);
	 PropagateExternal (inst, handleWE, l == LAYER_BY);

	 prev = inst;
      } /* for inst */
   }    /* for layer */
  
   this->RelaxBeamLimit();

   ++_dec->stats.nFrames;

   _dec->outPCache->cacheHit = _dec->outPCache->cacheMiss = 0;
   _dec->lmCache->transHit = _dec->lmCache->transMiss = 0;
   _dec->lmCache->laHit = _dec->lmCache->laMiss = 0;

   if (_dec->nPhone > 0)
      CalcPhonePost ();
}


/*  CC-mode style info for emacs
 Local Variables:
 c-file-style: "htk"
 End:
*/
