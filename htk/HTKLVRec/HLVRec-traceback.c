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
/*         File: HLVRec-traceback.c Traceback for              */
/*                                  HTK LV Decoder             */
/* ----------------------------------------------------------- */

/* Print Path
 */
void Decoder::PrintPath (WordendHyp *we)
{
   Pron pron;

   for (; we; we = we->prev) {
      pron = _dec->net->pronlist[we->pron];
      if ((we->user & 3) == 1)
         pron = pron->next;             /* sp */
      else if ((we->user & 3) == 2)
         pron = pron->next->next;       /* sil */
      
      printf ("%s%d (%d %.3f %.3f)  ", pron->word->wordName->name, pron->pnum, 
              we->frame, we->score, we->lm);
   }
}

/* PrintTok
*/
void Decoder::PrintTok(Token *tok)
{
   printf ("lmState %p score %f path %p ", tok->lmState, tok->score, tok->path);
   PrintPath (tok->path);
   printf ("\n");
}

/* PrintRelTok
*/
void Decoder::PrintRelTok(RelToken *tok)
{
   printf ("lmState %p delta %f path %p ", tok->lmState, tok->delta, tok->path);
   PrintPath (tok->path);
   printf ("\n");
}

/* PrintTokSet
 */
void Decoder::PrintTokSet (TokenSet *ts)
{
   int i;

   printf ("n = %d  score = %f id = %d\n", ts->n, ts->score, ts->id);

   for (i = 0; i < ts->n; ++i) {
      printf ("  %d ", i);
      PrintRelTok (&ts->relTok[i]);
   }
}

/* BestTokSet

     returns best token set in network
*/

TokenSet* Decoder::BestTokSet ()
{
   LexNodeInst *bestInst;
   TokenSet *tsi, *ts;
   LogFloat best;
   int i, N;

   bestInst = _dec->bestInst;
   switch (bestInst->node->type) {
   case LN_MODEL:
      N = bestInst->node->data.hmm->numStates;
      break;
   case LN_CON:
   case LN_WORDEND:
      N = 1;
      break;
   default:
      abort ();
      break;
   }
   
   ts = NULL;
   best = LZERO;
   for (i = 0; i < N; ++i) {
      tsi = &bestInst->ts[i];
      if (tsi->n > 0 && tsi->score > best)
         ts = tsi;
   }
   return (ts);
}

/* TraceBack

     Finds best token in end state and returns path.

*/
Transcription* Decoder::TraceBack(MemHeap *heap)
{
   Transcription *trans;
   LabList *ll;
   LLink lab, nextlab;
   WordendHyp *weHyp;
   TokenSet *ts;
   RelToken *bestTok;
   LogFloat prevScore, score;
   RelTokScore bestDelta;
   Pron pron;
   int i;
   HTime start;

   if (_dec->net->end->inst && _dec->net->end->inst->ts->n > 0)
      ts = _dec->net->end->inst->ts;
   else {
      HError (-9999, "no token survived to sent end!");

      ts = BestTokSet ();
      if (!ts) {        /* return empty transcription */
         HError (-9999, "best inst is dead as well!");
         trans = CreateTranscription (heap);
         ll = CreateLabelList (heap, 0);
         AddLabelList (ll, trans);
         
         return trans;
      }
   }

   bestDelta = LZERO;
   for (i = 0; i < ts->n; ++i)
      if (ts->relTok[i].delta > bestDelta) {
         bestTok = &ts->relTok[i];
         bestDelta = bestTok->delta;
      }
   assert (bestDelta <= 0.1);   /* 0.1 for accuracy reasons */
   
   if (trace & T_PROP)
      PrintTokSet (ts);

   if (trace & T_TOP) {
      printf ("best score %f\n", ts->score);
      PrintRelTok (bestTok);
   }

   trans = CreateTranscription (heap);
   ll = CreateLabelList (heap, 0);

   /* going backwards from </s> to <s> */
   for (weHyp = bestTok->path; weHyp; weHyp = weHyp->prev) {
      lab = CreateLabel (heap, ll->maxAuxLab);
      pron = _dec->net->pronlist[weHyp->pron];
      if ((weHyp->user & 3) == 1)
         pron = pron->next;             /* sp */
      else if ((weHyp->user & 3) == 2)
         pron = pron->next->next;       /* sil */

      lab->labid = pron->outSym;
      lab->score = weHyp->score;
      lab->start = 0.0;
      lab->end = weHyp->frame * _dec->frameDur * 1.0e7;
      lab->succ = ll->head->succ;
      lab->pred = ll->head;
      lab->succ->pred = lab->pred->succ = lab;
   }

   start = 0.0;
   prevScore = 0.0;
   for (lab = ll->head->succ; lab != ll->tail; lab = lab->succ) {
      lab->start = start;
      start = lab->end;
      score = lab->score - prevScore;
      prevScore = lab->score;
      lab->score = score;
   }

   for (lab = ll->head->succ; lab != ll->tail; lab = nextlab) {
      nextlab = lab->succ;
      if (!lab->labid)          /* delete words with [] outSym */
         DeleteLabel (lab);
   }

   AddLabelList (ll, trans);
   
   return trans;
}

/* LatTraceBackCount

     recursively assign numbers to wordendHyps (lattice nodes) and at the 
     same time count weHyps + altweHyps (lattice links)
*/
void Decoder::LatTraceBackCount (WordendHyp *path, int *nnodes, int *nlinks)
{
   AltWordendHyp *alt;

   if (!path)
      return;

   /* the pronvar is encoded in the user field: 
      0: -  1:sp  2: sil  */

   if (path->user < 4) {      /* not seen yet */

      /*      path->score *= -1.0; */

      ++(*nlinks);
      LatTraceBackCount (path->prev, nnodes, nlinks);
      for (alt = path->alt; alt; alt = alt->next) {
         ++(*nlinks);
         LatTraceBackCount (alt->prev, nnodes, nlinks);
      }
      ++(*nnodes);
      path->user += *nnodes * 4;         /* preserve pronvar bits */

   }
}

/* Paths2Lat

     recursively create nodes and arcs for weHyp end and pre_decInstessors
*/
void Decoder::Paths2Lat (Lattice *lat, WordendHyp *path, int *na)
{
   int s, n;
   AltWordendHyp *alt;
   TokScore prevScore;
   LNode *ln;
   LArc *la;
   Pron pron;

   if (!path)
      return;

   n = (int) (path->user / 4);  /* current node (end node of arcs) */

   if (!lat->lnodes[n].hook) {      /* not seen yet */
      ln = &lat->lnodes[n];
      ln->hook = (Ptr) path;

      ln->n = n;
      ln->time = path->frame*_dec->frameDur;   /* fix frame duration! */
      pron = _dec->net->pronlist[path->pron];

      if ((path->user & 3) == 1)
         pron = pron->next;             /* sp */
      else if ((path->user & 3) == 2)
         pron = pron->next->next;       /* sil */
      else 
         assert ((path->user & 3) == 0);
               

      ln->word = pron->word;
      ln->v = pron->pnum;

      if (trace & T_LAT)
         printf ("I=%d t=%.2f W=%d\n", n, path->frame*_dec->frameDur, path->pron);

      
      la = &lat->larcs[*na];
      ++(*na);

      s = path->prev ? (int) (path->prev->user / 4) : 0;
      prevScore = path->prev ? path->prev->score : 0.0;

      la->start = &lat->lnodes[s];
      la->end = &lat->lnodes[n];
      /* add to linked lists  foll/pred */
      la->farc = la->start->foll;
      la->start->foll = la;
      la->parc = la->end->pred;
      la->end->pred = la;

      la->prlike = pron->prob;
      la->aclike = path->score - prevScore - path->lm 
         - la->prlike * _dec->pronScale;
      la->lmlike = (path->lm - _dec->insPen) / _dec->lmScale;

#ifdef MODALIGN
      if (_dec->modAlign) {
         int startFrame;

         startFrame = path->prev ? path->prev->frame : 0;

         la->lAlign = LAlignFromModpath (lat->heap, path->modpath,
                                         startFrame, &la->nAlign);
      }
#endif

      if (trace & T_LAT)
         printf ("J=%d S=%d E=%d a=%f l=%f\n", *na, s, n, la->aclike, la->lmlike);

      Paths2Lat (lat, path->prev, na);

      /* alternatives */
      for (alt = path->alt; alt; alt = alt->next) {

         la = &lat->larcs[*na];
         ++(*na);

         s = alt->prev ? (int) (alt->prev->user / 4) : 0;
         prevScore = alt->prev ? alt->prev->score : 0.0;

         la->start = &lat->lnodes[s];
         la->end = &lat->lnodes[n];
         /* add to linked lists  foll/pred */
         la->farc = la->start->foll;
         la->start->foll = la;
         la->parc = la->end->pred;
         la->end->pred = la;

         la->prlike = pron->prob;
         la->aclike = (path->score + alt->score) - prevScore - alt->lm - 
            la->prlike * _dec->pronScale;
         la->lmlike = (alt->lm - _dec->insPen) / _dec->lmScale;
         
#ifdef MODALIGN
         if (_dec->modAlign) {
            int startFrame;
            
            startFrame = alt->prev ? alt->prev->frame : 0;
            
            la->lAlign = LAlignFromAltModpath (lat->heap, alt->modpath, path->modpath,
                                               startFrame, &la->nAlign);
         }
#endif
         if (trace & T_LAT)
            printf ("J=%d S=%d E=%d a=%f l=%f\n", *na, s, n, la->aclike, la->lmlike);
         
         Paths2Lat (lat, alt->prev, na);
      }
   }
}

/* LatTraceBack

     produce Lattice from the wordEnd hypotheses recoded in _dec
*/
Lattice* Decoder::LatTraceBack (MemHeap *heap)
{
   Lattice *lat;
   int i, nnodes = 0, nlinks = 0;
   WordendHyp *sentEndWE;

   if (!_dec->net->end->inst)
      HError (-9999, "LatTraceBack: end node not active");
   else
      printf ("found %d tokens in end state\n", _dec->net->end->inst->ts->n);

   if (buildLatSE && _dec->net->end->inst && _dec->net->end->inst->ts->n == 1)
      sentEndWE = _dec->net->end->inst->ts->relTok[0].path;
   else {
      if (buildLatSE)
         HError (-9999, "no tokens in sentend -- falling back to BUILDLATSENTEND = F");
      sentEndWE = BuildLattice ();
   }

   if (!sentEndWE) {
      HError (-9999, "LatTraceBack: no active sil wordend nodes");
      if (forceLatOut) {
         HError (-9999, "LatTraceBack: forcing lattice output");
#ifdef MODALIGN
         if (_dec->modAlign) 
/*             HError (-9999, "LatTraceBack: forced lattice output not supported with model-alignment"); */
            sentEndWE = BuildForceLat ();
         else 
#endif
            sentEndWE = BuildForceLat ();
      }
   }
   if (!sentEndWE)
      return NULL;
   
   /* recursively number weHyps (nodes), count weHyp + altweHyp (links) */
   LatTraceBackCount (sentEndWE, &nnodes, &nlinks);

   ++nnodes;    /* !NULL lattice start node */
   printf ("nnodes %d nlinks %d\n", nnodes, nlinks);

   /*# create lattice */
   lat = NewLattice (heap, nnodes, nlinks);

   /* #### fill in info (e.g. lmscale, inspen, models) */
   lat->voc = _dec->net->voc;
   lat->utterance = _dec->utterFN;
   lat->vocab = _dec->net->vocabFN;
   lat->hmms = _dec->hset->mmfNames ? _dec->hset->mmfNames->fName : NULL;
   lat->net = _dec->lm->name;
   lat->lmscale = _dec->lmScale;
   lat->wdpenalty = _dec->insPen;
   lat->prscale = _dec->pronScale;
   lat->framedur = 1.0;
      
   for (i = 0; i < nnodes; ++i)
      lat->lnodes[i].hook = NULL;

   {
      int na;
      na = 0;
      /* create lattice nodes & arcs */
      Paths2Lat (lat, sentEndWE, &na);
   }
   
#ifdef MODALIGN
   if (_dec->modAlign)
      CheckLAlign (lat);
#endif
   return lat;
}



/************      model-level traceback */

#ifdef MODALIGN
LAlign* Decoder::LAlignFromModpath (MemHeap *heap, ModendHyp *modpath, int wordStart, short *nLAlign)
{
   ModendHyp *m;
   LAlign *lalign;
   MLink ml;
   int n;
   int startFrame = 0;
   
   n = 0;
   for (m = modpath; m; m = m->prev)
      if (m->ln->type == LN_MODEL)
         ++n;
   
   lalign = (LAlign*) New (heap, n * sizeof(LAlign));
   *nLAlign = n;
   
   for (m = modpath; m; m = m->prev) {
      if (m->ln->type == LN_MODEL) {
         startFrame = (m->prev ? m->prev->frame : wordStart);
         ml = FindMacroStruct (_dec->hset, 'h', (Ptr) m->ln->data.hmm);
         if (!ml)
            HError (9999, "LAlignFromModpath: model not found!");

         assert (m->frame >= startFrame);
         --n;
         lalign[n].state = -1;
         lalign[n].like = 0.0;
         lalign[n].dur = (m->frame - startFrame) * _dec->frameDur;
         lalign[n].label = ml->id;
      }
   }
   assert (n == 0);

   return (lalign);
}

LAlign *Decoder::LAlignFromAltModpath (MemHeap *heap, ModendHyp *modpath, ModendHyp *mainModpath, int wordStart, short *nLAlign)
{
   ModendHyp *m, *nextM;
   LAlign *lalign;
   MLink ml;
   int n;
   int startFrame = 0;
   
   /* check for WE in main modpath */
   n = 0;
   for (m = mainModpath; m; m = m->prev)
      if (m->ln->type == LN_WORDEND)
         break;
      else {
         if (m->ln->type == LN_MODEL)
            ++n;
      }
   
   /* if there are no WE models in main modpath then we are looking at
      the </s> link and should just call the normal LAlignFromModpath() */
   if (!m)
      return LAlignFromModpath (heap, modpath, wordStart, nLAlign);

   /* take the n first model entries from the main modpaths (upto the WE node)
      and then switch to the alt modpath */

   for (m = modpath; m; m = m->prev)
      if (m->ln->type == LN_MODEL)
         ++n;
   
   lalign = (LAlign*) New (heap, n * sizeof(LAlign));
   *nLAlign = n;
   
   for (m = mainModpath; m; m = nextM) {
      nextM = m->prev;
      if (m->ln->type == LN_MODEL) {
         startFrame = (m->prev ? m->prev->frame : wordStart);
         ml = FindMacroStruct (_dec->hset, 'h', (Ptr) m->ln->data.hmm);
         if (!ml)
            HError (9999, "LAlignFromModpath: model not found!");

         assert (m->frame >= startFrame);
         --n;
         assert (n >= 0);
         lalign[n].state = -1;
         lalign[n].like = 0.0;
         lalign[n].dur = (m->frame - startFrame) * _dec->frameDur;
         lalign[n].label = ml->id;
      }
      else if (m->ln->type == LN_WORDEND)
         nextM = modpath;
   }
   assert (n == 0);

   return (lalign);
}

void Decoder::PrintModPath (ModendHyp *m)
{
   MLink ml;
   char *s, *t;

   for (; m; m = m->prev) {
      s = "?";
      switch (m->ln->type) {
      case LN_WORDEND:
         t = "WE";
         s = _dec->net->pronlist[m->ln->data.pron]->outSym->name;
         break;
      case LN_CON:
         t = "CON";
         s = "";
         break;
      case LN_MODEL:
         t = "MOD";
         ml = FindMacroStruct (_dec->hset, 'h', (Ptr) m->ln->data.hmm);
         if (ml)
            s = ml->id->name;
      }
      printf ("(%d %s %s) ", m->frame, t, s);
   }
   printf ("\n");
}

/* Faking sentence end arc model alignment */
void FakeSEModelAlign(Lattice *lat, LArc *la)
{  
   la->nAlign = 1;   
      
   la->lAlign = (LAlign*) New (lat->heap, sizeof(LAlign));
   la->lAlign->state = -1;
   la->lAlign->dur = la->end->time - la->start->time;
   la->lAlign->label = GetLabId("sil", FALSE);
}

void Decoder::CheckLAlign (Lattice *lat)
{
   int i, j;
   LArc *la;
   float dur, laDur;
   Pron pron;

   for (i = 0, la = lat->larcs; i < lat->na; ++i, ++la) {
      if (la->nAlign == 0 || !la->lAlign) {
         if (forceLatOut) {
            /* Faking sentence end arc model alignment */
            FakeSEModelAlign(lat, la);
         }
         else {
            HError (9999, "CheckLAlign: empty model alignment for arc %d", i);
         }
      }

      for (pron = la->end->word->pron; pron; pron = pron->next)
         if (pron->pnum == la->end->v)
            break;
      assert (pron);

      laDur = (la->end->time - la->start->time);
      dur = 0.0;
      for (j = 0; j < la->nAlign; ++j) {
         dur += la->lAlign[j].dur;

      }

      if (fabs (dur - laDur) > _dec->frameDur/2)
         printf ("CheckLAlign: MODALIGN Sanity check failed! %d laDur %.2f  dur %.2f\n", i, laDur, dur);
   }
}
#endif




/* FakeSEpath

     helper functions for BuildLattice and BuildForceLat.
     takens token and add LM transition to </s>

*/

AltWordendHyp* Decoder::FakeSEpath (RelToken *tok, bool useLM)
{
   AltWordendHyp *alt = NULL;
   PronId endPronId;
   LMState dest;
   LMTokScore lmScore;
   
   endPronId = _dec->net->end->data.pron;

   if (useLM)
      lmScore = LMCacheTransProb (_dec->lm, tok->lmState, endPronId, &dest);
   else
      lmScore = 0.0;
   if (lmScore > LSMALL) {  /* transition for END state possible? */
      assert (!useLM || dest == (Ptr) 0xfffffffe);
      lmScore += _dec->insPen;

      alt = (AltWordendHyp *) New (&_dec->altweHypHeap, sizeof (AltWordendHyp));
      alt->next = NULL;

      if (!_dec->fastlmla) {
         assert (lmScore <= tok->lmscore + 0.1); /* might not be true for more aggressive LMLA? */
      }
      /* temporarily store full score in altWEhyp */
      alt->score = tok->delta + (lmScore - tok->lmscore);
      alt->lm = lmScore;
      alt->prev = tok->path;
   }
   
   return alt;
}

/* AltPathList2Path

     Create full WordendHyp with alternatives from list of AltWordendHyps

*/
WordendHyp *AltPathList2Path (DecoderInst *_dec, AltWordendHyp *alt, PronId pron)
{
   WordendHyp *path;
   AltWordendHyp *bestAlt, *a;
   TokScore bestAltScore = LZERO;
   AltWordendHyp **pAlt;
   int i;

   /* find best */
   for (a = alt; a; a = a->next)
      if (a->score > bestAltScore) {
         bestAltScore = a->score;
         bestAlt = a;
      }
   
   /* create full WordendHyp for best */
   path = (WordendHyp *) New (&_dec->weHypHeap, sizeof (WordendHyp));
   path->prev = bestAlt->prev;
   path->pron = pron;
   path->frame = _dec->frame;
   path->score = bestAlt->score;
   path->lm = bestAlt->lm;
   path->user = 0;

   i = 0;
   pAlt = &path->alt;
   for ( ; alt; alt = alt->next) {
      if (alt != bestAlt) {
         ++i;
         *pAlt = alt;
         pAlt = &alt->next;
         alt->score = alt->score - path->score;
      }
   }
   *pAlt = NULL;

   printf ("found %d arcs\n", i);
   return path;
}


/* BuildLattice

     construct WordendHyp structure at the end of a sentence from all
     the tokensets in the final state of the sil Nodes in the SIL layer.
*/
WordendHyp* Decoder::BuildLattice ()
{
   int N, i;
   LexNodeInst *inst;
   LexNode *ln;
   TokenSet *ts;
   HLink hmmSP;
   AltWordendHyp *alt, *altPrev;
   WordendHyp *path;
   RelToken *tok;
#ifdef MODALIGN
   ModendHyp *silModend = NULL;
#endif

   alt = altPrev = NULL;
   hmmSP = _dec->net->hmmSP;
   for (inst = _dec->instsLayer[LAYER_SIL]; inst; inst = inst->next) {
      ln = inst->node;
      if (ln->type == LN_MODEL && ln->data.hmm != hmmSP) {
         N = ln->data.hmm->numStates;
         ts = &inst->ts[N-1];

         for (i = 0; i < ts->n; ++i) {
            tok = &ts->relTok[i];
            /* we have to update the path's score & frame, since 
               UpdateWordEndHyp() never got called on this path,
               A side effect is that the </s> link will have 
               aclike=0.0 and 1 frame length*/
            tok->path->score = ts->score + tok->delta;
            tok->path->frame = _dec->frame - 1;
            
#ifdef MODALIGN
            /* skip the final (MOD sil) modpath entry.
               If the token is outside the beam we will not have added a (MOD SIL) entry
               in PropagateExternal(). Nevertheless we should use this token, as it might
               slip into the LATPRUNEBEAM due to the final </s> LM transition applied to
               all tokens. 
            */
            if (_dec->modAlign) {
               if (tok->modpath->ln == ln)
                  tok->modpath = tok->modpath->prev;
               tok->path->frame = tok->modpath->frame;
               tok->path->modpath = tok->modpath;
               
               if (!silModend) {
                  silModend = (ModendHyp*) New (&_dec->modendHypHeap, sizeof (ModendHyp));
                  silModend->frame = _dec->frame;
                  silModend->ln = ln;   /* dodgy, but we just need ln with 'sil' model... */
                  silModend->prev = NULL;
               }
            }
#endif
            
            alt = FakeSEpath (tok, TRUE);
            if (alt) {
               alt->score += ts->score;
               alt->next = altPrev;
               altPrev = alt;
#ifdef MODALIGN
               alt->modpath = silModend;
#endif
            }
         } /* for tok */
      }
   } /* for inst */

   if (!alt)   /* no token in sil models at all */
      return NULL;

   path = AltPathList2Path (_dec, alt, _dec->net->end->data.pron);
#ifdef MODALIGN
   path->modpath = silModend;
#endif

   return path;
}

AltWordendHyp* Decoder::BuildLatAltList (TokenSet *ts, bool useLM)
{
   AltWordendHyp *alt, *altPrev;
   RelToken *tok;
   int i;
#ifdef MODALIGN
   ModendHyp *silModend = NULL;
#endif

   alt = altPrev = NULL;
   for (i = 0; i < ts->n; ++i) {
      tok = &ts->relTok[i];

      alt = FakeSEpath (tok, useLM);
      if (alt) {
         alt->score += ts->score;
         alt->next = altPrev;
         altPrev = alt;
#ifdef MODALIGN
         alt->modpath = silModend;
#endif
      }
   }
   return alt;
}


WordendHyp *Decoder::BuildForceLat ()
{
   TokenSet *ts;
   WordendHyp *path;
   AltWordendHyp *alt;
   RelToken *tok;
   int i;
#ifdef MODALIGN
   ModendHyp *silModend = NULL;
#endif

   ts = BestTokSet ();
   for (i = 0; i < ts->n; ++i) {
      tok = &ts->relTok[i];
      /* we have to update the path's score & frame, since 
         UpdateWordEndHyp() never got called on this path,
         A side effect is that the </s> link will have 
         aclike=0.0 and 1 frame length*/
      if (tok->path) {
         tok->path->score = ts->score + tok->delta;
         tok->path->frame = _dec->frame - 1;
      }
   }

   alt = BuildLatAltList (ts, TRUE);

   
   if (!alt) {  /* no valid LM transitions, try without */
      HError (-9999, "BuildForceLat: no tokens survived with valid LM transitions, inserting LM 0.0 arcs.");
      alt = BuildLatAltList (ts, FALSE);
   }

   if (!alt) {   /* how can this happen? */
      HError (-9999, "BuildForceLat: unable to force building lattice, giving up. THIS SHOULDN'T HAPPEN!");
      return NULL;
   }

   path = AltPathList2Path (_dec, alt, _dec->net->end->data.pron);
#ifdef MODALIGN
   path->modpath = silModend;
#endif
   return path;
}
