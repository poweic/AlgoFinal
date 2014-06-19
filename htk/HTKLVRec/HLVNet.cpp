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
/*         File: HLVNet.c Network handling for HTK LV Decoder  */
/* ----------------------------------------------------------- */

char *hlvnet_version = "!HVER!HLVNet:   3.4.1 [GE 12/03/09]";
char *hlvnet_vc_id = "$Id: HLVNet.c,v 1.1.1.1 2006/10/11 09:54:55 jal58 Exp $";

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HWave.h"
#include "HLabel.h"
#include "HAudio.h"
#include "HParm.h"
#include "HDict.h"
#include "HModel.h"
/* all the above are necessary just to include HModel.h */
/* #### sort out  include file dependencies! */

#include "config.h"

#include "HLVLM.h"      /* for LMId */
#include "HLVNet.h"

#include <assert.h>

#define LIST_BLOCKSIZE 70

#undef DEBUG_LABEL_NET

/* ----------------------------- Trace Flags ------------------------- */

#define T_TOP 0001         /* top level Trace  */
#define T_NET 0002         /* network summary*/
#define T_NETCON 0004    /* network construction  */

static int trace=0;
static ConfParam *cParm[MAXGLOBS];      /* config parameters */
static int nParm = 0;

/* --------------------------- Initialisation ---------------------- */

/* EXPORT->InitLVNet: register module & set configuration parameters */
void InitLVNet(void)
{
   int i;
   
   Register(hlvnet_version,hlvnet_vc_id);
   nParm = GetConfig("HLVNET", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
   }
}

/* --------------------------- the real code  ---------------------- */

TLexNode::TLexNode(TLexNet* net, int layerId, HLink hmm): type(LN_MODEL) {
    this->init(net, layerId);
    this->hmm = hmm;
}

TLexNode::TLexNode(TLexNet* net, int layerId, Pron pron): type(LN_WORDEND) {
    this->init(net, layerId);
    this->pron = pron;
}

TLexNode::TLexNode(TLexNet* net, int layerId, LabId lc, LabId rc): type(LN_CON) {
    this->init(net, layerId);
    this->lc = lc;
    this->rc = rc;
}

void TLexNode::init(TLexNet* net, int layer_id) {
  nlinks = 0;
  link = NULL;
  ln = NULL;
  loWE = hiWE = lmlaIdx = 0;
  layerId = layer_id;

  chain = net->root;
  net->root = this;
  ++net->nNodesLayer[layerId];

  hmm = NULL;
  pron = NULL;
  lc = rc = NULL;
}

void* TLexNode::operator new (size_t, MemHeap* heap) {
  assert (heap);
  return New (heap, sizeof (TLexNode));
}

/* from HNet.c */
/* Binary search to find elem in n element array */
/* Compare to BSearch, this function add element to array if not found */
static int BAddSearch (Ptr elem, int &np, Ptr **ap)
{
   Ptr *array;
   int l,u,c;

   // Binary search
   array=*ap;
   l=1;u=np;
   while(l<=u) {
      c=(l+u)/2;
      if (array[c]==elem) return(c);
      else if (array[c]<elem) l=c+1; 
      else u=c-1;
   }

   // Add elem to array if NOT FOUND
   if (((np+1) % LIST_BLOCKSIZE)==0) {
      Ptr *newId;

      newId=(Ptr *) New(&gcheap,sizeof(Ptr)*(np+1+LIST_BLOCKSIZE));
      for (c=1;c<=np;c++)
         newId[c]=array[c];
      Dispose(&gcheap,array);
      *ap=array=newId;
   }
   for (c=1;c<=np;c++)
      if (elem<array[c]) break;

   for (u=np++;u>=c;u--)
      array[u+1]=array[u];

   array[c]=elem;
   return(c);
}

#define HASH1(p1,s)  (((int)(p1))%(s))
#define HASH2(p1,p2,s)  (( HASH1(p1,s) + HASH1(p2,s)) % (s))

TLexConNode* TLexNet::FindAddTLCN (LayerId layerId, int *n,
    TLexConNode *lcnHashTab[], LabId lc, LabId rc)
{
   unsigned int hash = HASH2(lc, rc, LEX_CON_HASH_SIZE);

   TLexConNode *lcn = lcnHashTab[hash];
   while (lcn && !((lcn->lc == lc) && (lcn->rc == rc)) )
      lcn = lcn->next;

   if (!lcn) {          /* not found */
     assert(heap);
      // if (!heap) return NULL;
      lcn = (TLexConNode *) new(heap) TLexNode( this, layerId, lc, rc);

      lcn->next = lcnHashTab[hash];
      lcnHashTab[hash] = lcn;
      ++(*n);
   }
   return lcn;
}

/* add TLexNode to hashtable
   only hmm is used for hashing and comparison, 
   i.e. no two node with same hmm, but different types are possible
*/
TLexNode* TLexNet::FindAddTLexNode (LayerId layerId, int *n,
    TLexNode *lnHashTab[], LexNodeType type , HLink hmm)
{
   unsigned int hash = HASH1 (hmm, LEX_MOD_HASH_SIZE);

   TLexNode *ln = lnHashTab[hash];
   while (ln && !((ln->hmm == hmm)) )
      ln = ln->next;

   /* not found */
   if (!ln) {
      if (!heap) return NULL;
      ln = new(heap) TLexNode( this, layerId, hmm);

      ln->next = lnHashTab[hash];
      lnHashTab[hash] = ln;
      ++(*n);
   }
   return ln;
}

HLink FindTriphone (HMMSet *hset, LabId a, LabId b, LabId c)
{
   const size_t BUFLEN = 100;
   char buf[BUFLEN];
   
   if (sprintf (buf, "%s-%s+%s", a->name, b->name, c->name) == -1) 
      HError (9999, "HLVNet: model names too long");
   
   LabId triLabId = GetLabId (buf, FALSE);
   if (!triLabId)
      HError (9999, "HLVNet: no model label for phone (%s)", buf);
   
   MLink triMLink = FindMacroName (hset, 'l', triLabId);
   if (!triMLink)
      HError (9999, "HLVNet: no model macro for phone (%s)", buf);
   
   return ((HLink) triMLink->structure);
}

/*  scan vocabulary pronunciations for phone sets A, B, AB, YZ
    sets A and B are small (~50) and stored in arays
    sets AB and YZ are larger (~750) and are stored in hash tables.
*/
void TLexNet::CollectPhoneStats ()
{
   LabId sil, z, p;

   for (int i = 0; i < VHASHSIZE; i++) {
     for (auto word = voc->wtab[i]; word ; word = word->next) {
       if (word->aux != (Ptr) 1 || (word->wordName == startId || word->wordName == endId) )
	 continue;
       
       // word->print();
       for (auto pron = word->pron; pron ; pron = pron->next) {
	 // pron->print();
	 const auto& N = pron->nphones;
	 const auto& phones = pron->phones;

	 if (N < 1)
	   HError (9999, "CollectPhoneStats: pronunciation of '%s' is empty", word->wordName->name);

	 if (pron->aux != (Ptr) 1)
	   continue;

	 add_phones (phones[ 0 ], 'A');
	 add_phones (phones[N-1], 'Z');

	 if (N == 1) {
	   add_phones (phones[0], 'P');
	   /*#### need to add ZP node in YZ list for single phone word */
	   /*#### here only collect phones used in single phone words */
	   continue;
	 }

	 /* add to AB and YZ hashes */
	 FindAddTLCN (LAYER_AB, &nlexAB, lexABhash, phones[ 0 ], phones[ 1 ]);
	 FindAddTLCN (LAYER_YZ, &nlexYZ, lexYZhash, phones[N-2], phones[N-1]);

       }
     }
   }
   
   /* add 'sil' to A and Z lists */
   sil = GetLabId ("sil", FALSE);
   if (!sil)
      HError (9999, "cannot find 'sil' model.");

   add_phones (sil, 'A');
   add_phones (sil, 'Z');

   /* for each phone P occuring in a single phone word 
        for each word end context Z
          add node ZP to SA and YZ layers
   */

   if (trace & T_NETCON)
      printf ("adding extra nodes for single phone words:\n");

   for (int i = 1; i <= nlexP; ++i) {
      if (trace & T_NETCON)
         printf ("P='%s'  \33[33m => \33[0m ", lexP[i]->name);
      p = lexP[i];
      for (int j = 1; j <= nlexZ; ++j) {
         if (trace & T_NETCON)
            printf ("z='%s' ", lexZ[j]->name);
         z = lexZ[j];

         FindAddTLCN (LAYER_SA, &nlexSA, lexSAhash, z, p);
         FindAddTLCN (LAYER_YZ, &nlexYZ, lexYZhash, z, p);
      }
      if (trace & T_NETCON)
         printf ("\n");
   }

   if (trace & T_NETCON) { 
      TLexNode *lcn;

      /* debug: print lexA, lexZ, lexAB & lexYZ*/
      printf ("\33[33mnlexA = %d   \33[0m", nlexA);
      for (int i = 1; i <= nlexA; ++i)
         printf ("%s ", lexA[i]->name);
      printf ("\n");
      
      printf ("\33[33mnlexZ = %d   \33[0m", nlexZ);
      for (int i = 1; i <= nlexZ; ++i)
         printf ("%s ", lexZ[i]->name);
      printf ("\n");
      
      printf ("\33[33mnlexAB = %d   \33[0m", nlexAB);
      for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
         for (lcn = lexABhash[i]; lcn; lcn = lcn->next)
            printf ("%s-%s  ", lcn->lc->name, lcn->rc->name);
      printf ("\n");
      
      printf ("\33[33mnlexYZ = %d   \33[0m", nlexYZ);
      for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
         for (lcn = lexYZhash[i]; lcn; lcn = lcn->next)
            printf ("%s-%s  ", lcn->lc->name, lcn->rc->name);
      printf ("\n");
   }

}

/* add link from TLexNode start to end, if it doesn't exist already 
*/
void AddLink (MemHeap *heap, TLexNode *start, TLexNode *end)
{
   TLexLink *ll;

   assert (start);
   assert (end);
   /* check if link exists already */
   /*#### GE: make check optional? */
   for (ll = start->link; ll; ll = ll->next)
      if (ll->end == end)
         break;

   if (!ll) {
      ll = (TLexLink *) New (heap, sizeof (TLexLink));
      ll->start = start;
      ll->end = end;
      
      ++start->nlinks;
      ll->next = start->link;
      start->link = ll;
   }
}

/* Find link from node ln to a node corresponding to model hmm
 */
TLexLink *FindHMMLink (TLexNode *ln, HLink hmm)
{
   TLexLink *ll;

   for (ll = ln->link; ll; ll = ll->next)
      if (ll->end->hmm == hmm)
         break;

   return ll;
}

/* Create initial phone (A) layer of z-a+b nodes 
*/
void TLexNet::CreateAnodes ()
{     
   LabId z, a, b;
   
   /* create word initial (A) nodes */
   /*   for each AB node */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i) 
      for (auto lcnAB = lexABhash[i]; lcnAB; lcnAB = lcnAB->next) {
         a = lcnAB->lc;
         b = lcnAB->rc;
         /* for each Z node */
         for (int j = 1; j <= nlexZ; ++j) {
            z = lexZ[j];

            /* create Node z-a+b */
            HLink hmm = FindTriphone (hset, z, a, b);
            auto lnZAB = FindAddTLexNode (LAYER_A, &nNodeA, nodeAhash, LN_MODEL, hmm);
            lnZAB->lc = NULL;

            /* connect to ZA node */
            auto lcnZA = FindAddTLCN (LAYER_SA, &nlexSA, lexSAhash, z, a);
            AddLink (heap, lcnZA, lnZAB);

            /* connect to AB node */
            AddLink (heap, lnZAB, lcnAB);
         }
      }

   if (trace & T_NETCON) {
      printf ("nlexSA = %d   ", nlexSA);
      for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
         for (auto lcn = lexSAhash[i]; lcn; lcn = lcn->next)
            printf ("%s-%s  ", lcn->lc->name, lcn->rc->name);
      printf ("\n");

      printf ("nNodeA = %d   ", nNodeA);
      for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
         for (auto lnZAB = nodeAhash[i]; lnZAB; lnZAB = lnZAB->next)
            printf ("zab_hmm %p nlinks %d\n", lnZAB->hmm, lnZAB->nlinks);
   }
}


/* Create final phone (Z) layer of y-z+a nodes 
*/
void TLexNet::CreateZnodes ()
{
   LabId y, z, a;

   /* create word final (Z) nodes */
   /*   for each YZ node */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i) 
      for (auto lcnYZ = lexYZhash[i]; lcnYZ; lcnYZ = lcnYZ->next) {
         y = lcnYZ->lc;
         z = lcnYZ->rc;
         /* for each A phone */
         for (int j = 1; j <= nlexA; ++j) {
            a = lexA[j];
            /* create Node y-z+a */
            HLink hmm = FindTriphone (hset, y, z, a);
            auto lnYZA = FindAddTLexNode (LAYER_Z, &nNodeZ, nodeZhash, LN_MODEL, hmm);
            lnYZA->lc = NULL;

            /* connect to YZ node */
            AddLink (heap, lcnYZ, lnYZA);

            /* connect to ZS node */
            auto lcnZA = FindAddTLCN (LAYER_ZS, &nlexZS, lexZShash, z, a);
            AddLink (heap, lnYZA, lcnZA);
         }
      }
   if (trace & T_NETCON) {
      printf ("nlexZS = %d   ", nlexZS);
      for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
         for (auto lcn = lexZShash[i]; lcn; lcn = lcn->next)
            printf ("%s-%s  ", lcn->lc->name, lcn->rc->name);
      printf ("\n");

      printf ("nNodeZ = %d   ", nNodeZ);
      for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
         for (auto lnYZA = nodeZhash[i]; lnYZA; lnYZA = lnYZA->next)
            printf ("yza_hmm %p nlinks %d\n", lnYZA->hmm, lnYZA->nlinks);
   }
}

/* get HLink from LabId
   #### GE: this should really be in HModel
 */
HLink FindHMM (HMMSet *hset, LabId id)
{
   assert (id);

   MLink ml = FindMacroName (hset, 'l', id);
   if (!ml)
      HError (9999, "cannot find model for label '%s'", id->name);

   assert (ml->structure);
   return ((HLink) ml->structure);
}

/* Create sil/sp model nodes between ZS and SA nodes 
 */
void TLexNet::CreateSILnodes ()
{
  /* find sil & sp models */
  LabId sil = GetLabId ("sil", FALSE);
  if (!sil) HError (9999, "cannot find 'sil' model.");
  HLink hmmSIL = FindHMM (hset, sil);

  LabId sp = GetLabId ("sp", FALSE);
  if (!sp) HError (9999, "cannot find 'sp' model.");
  HLink hmmSP = FindHMM (hset, sp);

  for (int i = 0; i < LEX_CON_HASH_SIZE; ++i) {
    for (auto lcnZS = lexZShash[i]; lcnZS; lcnZS = lcnZS->next) {
      LabId z = lcnZS->lc;
      LabId s = lcnZS->rc;
      /*         printf ("ZS node %p %s-%s i %d\n", lcnZS, z->name, s->name, i); */
      auto ln = new(heap) TLexNode( this, LAYER_SIL, (s==sil) ? hmmSIL : hmmSP);

      ln->next = nodeSIL;
      nodeSIL = ln;
      ++nNodeSIL;

      AddLink (heap, lcnZS, ln);
      if (s == sil) {        /* sil node */
	ln->lc = sil;

	/* connect sil node to all sil-A nodes */
	for (int j = 1; j <= nlexA; ++j) {
	  if (lexA[j] == sil)
	    continue;

	  /* list of A contxts includes sil! */
	  auto lcnSA = FindAddTLCN (LAYER_SA, &nlexSA, lexSAhash, s, lexA[j]);
	  /*  printf ("  sil SA node %s-%s\n", lcnSA->lc->name, lcnSA->rc->name); */
	  AddLink (heap, ln, lcnSA);
	}
      }
      else {         /* sp node */
	ln->lc = sp;

	/* connect sp node to corresponding SA node (SA==ZS) */
	auto lcnSA = FindAddTLCN (LAYER_SA, &nlexSA, lexSAhash, z, s);
	/*  printf ("  sp  SA node %s-%s\n", lcnSA->lc->name, lcnSA->rc->name); */
	AddLink (heap, ln, lcnSA);
      }

    }
  }
}

void TLexNet::add_phones(LabId elem, char type) {
  switch (type) {
    case 'A':
      BAddSearch((Ptr) elem, nlexA, ((Ptr **) &lexA));
      break;
    case 'Z':
      BAddSearch((Ptr) elem, nlexZ, ((Ptr **) &lexZ));
      break;
    case 'P':
      BAddSearch((Ptr) elem, nlexP, ((Ptr **) &lexP));
      break;
  }
}

/* Handle1PhonePron

     create network structure for single phone words
     The actual phones are put in the Z layer.
     duplicate wordend for each preceeding context Z

     ZP ----> WE ---> ZP ----> z-p+  ---> PS
     SA               YZ       y-z+a ---> ZS    

*/
void TLexNet::Handle1PhonePron (Pron pron)
{
   assert (pron->nphones == 1);

   if (trace & T_NETCON)
      printf ("handle single phone pron for '%s'\n", pron->word->wordName->name);

   /* the single phone */
   LabId p = pron->phones[0];

   PronId pronid = ++this->nPronIds;
   int lmlaIdx = ++this->lmlaCount;

   pron->aux = (Ptr) (int) pronid;

   /* for each word end phone Z
        create a WE node and link to appropriate nodes in SA & YZ layers
        zp_sa ---> WE ---> zp_yz
   */
   for (int i = 1; i <= this->nlexZ; ++i) {
      if (trace & T_NETCON)
         printf ("Z node '%s' P node '%s'\n", this->lexZ[i]->name, p->name);
      LabId z = this->lexZ[i];

      auto zp_sa = FindAddTLCN (LAYER_SA, &this->nlexSA, this->lexSAhash, z, p);
      auto zp_yz = FindAddTLCN (LAYER_YZ, &this->nlexYZ, this->lexYZhash, z, p);

      /* create word end node */
      auto we = new(heap) TLexNode( this, LAYER_WE, pron);
      
      we->loWE = we->hiWE = pronid;
      we->lmlaIdx = lmlaIdx;

      we->lc = pron->word->wordName;
      we->next = this->nodeBY; /*#### link into BY list? */
      this->nodeBY = we;
      ++this->nNodeBY;

      AddLink (heap, zp_sa, we);
      AddLink (heap, we, zp_yz);
   }

}

/* CreateBYnodes

     Create the forward tree from phone B to phone Y, 
     linking the appropriate AB and YZ nodes.
*/
void TLexNet::CreateBYnodes () 
{
   int nshared = 0;

   /* Create tree B -- Y */

   /* for each pron with 2 or more phones */
   for (int i = 0; i < VHASHSIZE; i++) {
     for (auto word = voc->wtab[i]; word ; word = word->next) {
       if ( word->aux != (Ptr) 1 || (word->wordName == startId || word->wordName == endId) )
	 continue;

       for (auto pron = word->pron; pron ; pron = pron->next) {
	 if (pron->aux != (Ptr) 1) continue;

	 const auto& N = pron->nphones;
	 const auto& phones = pron->phones;

	 if (N < 2) {
	   Handle1PhonePron (pron);
	   continue;
	 }

	 /* find AB node (see origina HLVnet.c find only (no add) )*/
	 TLexNode* prevln = FindAddTLCN (LAYER_AB, &nlexAB, lexABhash, phones[0], phones[1]);

	 /* add models for phones B -- Y */
	 for (int p = 1; p < N - 1; ++p) {
	   HLink hmm = FindTriphone (hset, phones[p-1], phones[p], phones[p+1]);

	   /* search in prevln's successors */
	   TLexLink *ll = FindHMMLink (prevln, hmm);
	   if (ll) {
	     prevln = ll->end;
	     ++nshared;
	   }
	   else {             /* model not found -> create a new one */
	     auto ln = new(heap) TLexNode( this, LAYER_BY, hmm);

	     ln->next = nodeBY;
	     nodeBY = ln;
	     ++nNodeBY;

	     AddLink (heap, prevln, ln);   /* guaranteed to be a new link! */

	     ln->lc = NULL;
	     prevln = ln;
	   }
	 }
	 /* create word end node */
	 auto ln = new(heap) TLexNode( this, LAYER_WE, pron);

	 ln->lc = pron->word->wordName;
	 ln->next = nodeBY;
	 nodeBY = ln;
	 ++nNodeBY;

	 AddLink (heap, prevln, ln);   /* guaranteed to be a NEW link! */
	 prevln = ln;

	 /* find YZ node (see origina HLVnet.c find only (no add) )*/
	 ln = FindAddTLCN (LAYER_YZ, &nlexYZ, lexYZhash, phones[N - 2], phones[N - 1]);

	 /* find LexNode and connect prevln to it */
	 AddLink (heap, prevln, ln);

       } // End of pron for-loop
       
     }
   }
   if (trace & T_NETCON)
      printf ("nodes shared in prefix tree: %d\n", nshared);
}


/* CreateBoundary

     create one model and one word end node for a given boundary (start or end) labid.
     return WE node with lnWE->next pointin to model node
*/
TLexNode* TLexNet::CreateBoundary (LabId labid, int modLayer, int weLayer)
{
   Word w = GetWord (this->voc, labid, FALSE);

   if (!w)
      HError (9999, "HLVNet: cannot find START/ENDWORD '%s'", labid->name);

   if (w->nprons != 1 || w->pron->nphones != 1)
      HError (9999, "HLVNet: only one pronuinciation with one model allowed for START/ENDWORD '%s'", labid->name);
   
   /* create model node */
   auto lnMod = new(heap) TLexNode( this, modLayer, FindHMM (this->hset, w->pron->phones[0]));

   lnMod->next = NULL;
   lnMod->lc = w->pron->phones[0];

   /* create word end node */
   auto lnWe = new(heap) TLexNode( this, weLayer, w->pron);

   lnWe->lc = w->pron->word->wordName;
   lnWe->next = lnMod;
   
   return lnWe;
}

/* CreateStartEnd 

     create model and word end nodes for STARTWORD and ENDWORD and connect them 
     to the appropriate places.

*/
void TLexNet::CreateStartEnd ()
{
   TLexNode *lnWe, *lnMod, *lnTime;

   /* STARTWORD */
   /* the two layers must be before SA; do NOT use ZS (spSkipLayer) */   
   lnWe = CreateBoundary (this->startId, LAYER_Z, LAYER_SIL); 
   lnMod = lnWe->next;
   lnWe->next = NULL;
   AddLink (heap, lnMod, lnWe);

   lnWe->loWE = lnWe->hiWE = ++this->nPronIds;
   lnWe->lmlaIdx = ++this->lmlaCount;
   lnWe->pron->aux = (Ptr) lnWe->loWE;

   /* connect start WE -> all matching SA nodes */
   /* loop over all SA nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexSAhash[i]; ln; ln = ln->next)
         if (ln->lc == lnMod->lc) {         /* matching left context? */
            AddLink (heap, lnWe, ln);
         }

   this->start = lnMod;

   /* ENDWORD */
   /*   put null node (lnTime) into LAYER_SA to get updated end time
        in UpdateWordEndHyp() really only for HVite compatible times.
        this node is also used to detect the path to the Sentend silence for 
        the pronprob handling in HLVRec!  */
   
   lnWe = CreateBoundary (this->endId, LAYER_A, LAYER_AB);  
   lnMod = lnWe->next;
   lnWe->next = NULL;
   AddLink (heap, lnMod, lnWe);
   
   lnWe->loWE = lnWe->hiWE = ++this->nPronIds;
   lnWe->lmlaIdx = ++this->lmlaCount;
   lnWe->pron->aux = (Ptr) lnWe->loWE;
   
   lnTime = new(heap) TLexNode( this, LAYER_SA, lnWe->lc, lnWe->lc);
   AddLink (heap, lnTime, lnMod);
   
   /* connect start all matching ZS nodes -> LN time  */
   /* loop over all ZS nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexZShash[i]; ln; ln = ln->next)
         if (ln->rc == lnMod->lc) {         /* matching right context? */
            AddLink (heap, ln, lnTime);
         }

   this->end = lnWe;

   /* for -/sp/sil dicts add an extra sp and sil model leading to SENT_END */
   if (this->silDict) {
      /* find sil & sp models */
      LabId sil = GetLabId ("sil", FALSE);
      if (!sil)
         HError (9999, "cannot find 'sil' model.");
      HLink hmmSIL = FindHMM (this->hset, sil);

      LabId sp = GetLabId ("sp", FALSE);
      if (!sp)
         HError (9999, "cannot find 'sp' model.");
      HLink hmmSP = FindHMM (this->hset, sp);
      
      this->lnSEsp = new(heap) TLexNode( this, LAYER_SIL, hmmSP);
      this->lnSEsil = new(heap) TLexNode( this, LAYER_SIL, hmmSIL);

      AddLink (heap, this->lnSEsp, lnTime);
      AddLink (heap, this->lnSEsil, lnTime);
   }

}


/*  TraverseTree

      number all wordend nodes below ln, starting with Id start and return
      highest Id. In every node store lo and hi Ids of WE nodes
      in subtree below it.

*/
int TraverseTree (TLexNode *ln, int start, int &lmlaCount)
{
   int curHi;

   /* WE? then bottom out, get new id, store it with ln and return */
   if (ln->type == LN_WORDEND) {
      ln->loWE = ln->hiWE = start;    /* store ID in ln->pron */
      ln->pron->aux = (Ptr) ln->loWE;

      // printf("LMLA %p %d %d  %d\n", ln, ln->loWE, ln->hiWE, lmlaCount);
      return start;
   }

   ln->loWE = start;
   curHi = start - 1;
   for (auto ll = ln->link; ll; ll = ll->next) {
      start = curHi + 1;
      curHi = TraverseTree (ll->end, start, lmlaCount);
   }
   ln->hiWE = curHi;
   /* ^^^ One could probably replace curHi by ln->hiWE in the above ^^^ */
   
   /* unique successors don't get a slot in the LMlaTree */
   if (ln->nlinks > 1) {
      for (auto ll = ln->link; ll; ll = ll->next)
         ll->end->lmlaIdx = ++lmlaCount; // Assign new lmlaIds to successor 
      // lmlaCount += ln->nlinks;
   }
   else {
      /* LA Tree Compression :
       * no LM handling required for the unique successor of this node */

      /* LM look ahead tree.
      
           Each node in the prefix part of the LexNet points to an entry
           in LMlaTree. An entry of 0 indicates that the lookahead
           information doesn't need to be updated. This occur if the set
           of reachable word ends is the same as for the predecessor node,
           e.g.:
                     3-->0
                    /
           1-->0-->0
                    \
                     2
      */

      ln->link->end->lmlaIdx = 0;
   }

   // printf("LMLA %p %d %d  %d\n", ln, ln->loWE, ln->hiWE, lmlaCount);

   return curHi;
}

/* AssignWEIds
 *
 * */
void TLexNet::AssignWEIds()
{
   /* assign PronIds to wordend nodes reachable from AB nodes.
      this excludes:
       - single phone prons (cf. Handle1PhonePron)
       - start/end words    (cf. CreateStartEnd)
   */

   int highest = this->nPronIds;
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i) {
      for (auto ln = this->lexABhash[i]; ln; ln = ln->next) {
         int start = highest + 1;

         /*          ln->loWE = start; */
         ln->lmlaIdx = ++this->lmlaCount;
         highest = TraverseTree (ln, start, this->lmlaCount);
         ln->hiWE = highest;


         if (trace & T_NETCON)
            printf ("AB %s-%s  loWE %d hiWE %d  lmlaIdx %d\n", ln->lc->name, ln->rc->name, 
                    ln->loWE, ln->hiWE, ln->lmlaIdx);
      }
   }

   /* one A node can have multiple successors. */
   /* create complex lmla nodes for them in CreateCompLMLA() */


   /*# handle single phone prons */

   /*# handle <s> and </s> */

   /*# sanity check for prons without id */

   printf ("lmla count %d\n", this->lmlaCount );
   printf ("nprons %d\n", highest);

   if (sizeof(PronId) < 4 || sizeof (LMId) < 4) {
      if (highest > (1UL << (8 * sizeof (PronId))) ||
          highest > (1UL << (8 * sizeof (LMId))))
         HError (9999, "AssignWEIds: too many pronunciations for PronId/LMId type. Recompile with type int");
   }
}


void CreateCompLMLA (MemHeap *heap, LMlaTree *laTree, TLexNet *tnet)
{
   /* first handle simple case and count complex cases */
   int nComp = 0;
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
      for (auto ln = tnet->nodeAhash[i]; ln; ln = ln->next) {
         if (!ln->link->next) {      /* unique successor node */
            /* copy lo&hi to A nodes from AB nodes */
            ln->lmlaIdx = ln->link->end->lmlaIdx;
            ln->loWE = ln->link->end->loWE;
            ln->hiWE = ln->link->end->hiWE;
         }
         else {
            ++nComp;
         }
      }

   /* alloc space for comp nodes */
   laTree->nCompNodes = nComp;

   laTree->compNode = (CompLMlaNode *) New (heap, nComp * sizeof (CompLMlaNode));

   /* fill in info */
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i) {
      for (auto ln = tnet->nodeAhash[i]; ln; ln = ln->next) {
         if (!ln->link->next) {      /* simple case */
         }
         else {
            /* construct complex LMLA node */
            int nfoll = 0;
            for (auto tll = ln->link; tll; tll = tll->next)
               ++nfoll;

            /* create CompLMLANode */
            ln->lmlaIdx = ++tnet->lmlaCount;

            auto compNode = &laTree->compNode[ln->lmlaIdx - laTree->nNodes];
            compNode->n = nfoll;
            compNode->lmlaIdx = (int*) New (heap, nfoll * sizeof (int));
            
	    int j = 0;
            for (auto tll = ln->link; tll; tll = tll->next, ++j)
               compNode->lmlaIdx[j] = tll->end->lmlaIdx;
         }
      }
   }

   assert (tnet->lmlaCount+1 - laTree->nNodes == laTree->nCompNodes); /* extra entry for lmlaIdx=0 */
}

/*

*/
void InitLMlaTree(LexNet *net, TLexNet *tnet)
{
   /* simple nodes */
   LMlaTree *laTree = (LMlaTree *) New (net->heap, sizeof(LMlaTree));
   net->laTree = laTree;
   laTree->nNodes = tnet->lmlaCount + 1;        /* extra entry for lmlaIdx=0 */
   laTree->node = (LMlaNode *) New (net->heap, laTree->nNodes * sizeof (LMlaNode));

   /* complex nodes */
   CreateCompLMLA (net->heap, laTree, tnet);
}

MemHeap LexNet::_netHeap("Net heap", MSTAK, 1, 0,100000, 800000);

/* Create the Lexicon Network based on the vocabulary and model set
*/
void LexNet::init(Vocab *voc, HMMSet *hset, 
    char *startWord, char *endWord, bool silDict) {

   TLexNet* tnet = new TLexNet(voc, hset, startWord, endWord, silDict);
   tnet->init();

   /* convert TLexNet to more compact LexNet */
   this->ConvertTLex2Lex(tnet);

   /* all tokens pass through SA directly before (i.e. with no time diff) the first 
      model of a new word. Update time and score in last weHyp of token in this layer */
   this->wordEndLayerId = LAYER_SA;

   this->startPron = (PronId) (int) GetWord (voc, tnet->startId, FALSE)->pron->aux;
   this->endPron = (PronId) (int) GetWord (voc, tnet->endId, FALSE)->pron->aux;

   /* add pronprobs in ZS nodes and (if S==sp) propagate - variant bypassing sp model. */
   this->spSkipLayer = LAYER_ZS;

   LabId spLab = GetLabId ("sp", FALSE);
   if (!spLab)
     HError (9999, "cannot find 'sp' model.");
   this->hmmSP = FindHMM (this->hset, spLab);

   this->silDict = silDict;

   /* get rid of temporary data structures */
   TLexNet::ResetHeap();
}

void* LexNet::operator new (size_t) {
  return New (&LexNet::_netHeap, sizeof (LexNet));
}

/* ConvertTLex2Lex

     convert the large, verbose temp lex structure into the compact LexNet structure
     nodes are ordered in layers
 */
void LexNet::ConvertTLex2Lex (TLexNet *tnet)
{
   LexNode *layerCur[NLAYERS];        /* pointer to next free LN in layer */

   if (trace &T_NET) {
      printf ("number of nodes:\n");
      printf ("  A model nodes:    %d\n", tnet->nNodeA);
      printf ("  AB context nodes: %d\n", tnet->nlexAB);
      printf ("  B-Y model nodes:  %d\n", tnet->nNodeBY);
      printf ("  YZ context nodes: %d\n", tnet->nlexYZ);
      printf ("  Z model nodes:    %d\n", tnet->nNodeZ);
      printf ("  ZS context nodes: %d\n", tnet->nlexZS);
      printf ("  SIL model nodes:  %d\n", tnet->nNodeSIL);
      printf ("  SA context nodes: %d\n", tnet->nlexSA);
      printf ("total: %d\n", tnet->nNodeA + tnet->nlexAB + tnet->nNodeBY + tnet->nlexYZ +
              tnet->nNodeZ + tnet->nlexZS + tnet->nlexSA + tnet->nNodeSIL);
   }

   int nn = 0;
   for (int i = 0; i < NLAYERS; ++i) {
      nn += tnet->nNodesLayer[i];
      if (trace & T_NET)
         printf ("layer %d contains %d nodes\n", i, tnet->nNodesLayer[i]);
   }

   if (trace & T_NET)
      printf ("total: %d nodes\n", nn);
   
   if (tnet->silDict)
      /* 7 special models: <s> Mod/WE,  </s> Time, </s> Mod/WE, </s> SP, </s> SIL */
      assert (nn == tnet->nNodeA + tnet->nlexAB + tnet->nNodeBY + tnet->nlexYZ +
              tnet->nNodeZ + tnet->nlexZS + tnet->nlexSA + tnet->nNodeSIL + 7);
   else
      /* 5 special models: <s> Mod/WE,  </s> Time, </s> Mod/WE */
      assert (nn == tnet->nNodeA + tnet->nlexAB + tnet->nNodeBY + tnet->nlexYZ +
              tnet->nNodeZ + tnet->nlexZS + tnet->nlexSA + tnet->nNodeSIL + 5);
   

   /* alloc and init LexNet */
   // LexNet *net = (LexNet *) New (heap, sizeof (LexNet));
   this->heap = &LexNet::_netHeap;
   this->voc = tnet->voc;
   this->hset = tnet->hset;
   this->nLayers = NLAYERS;      /* keep HLVRec general without ugly constants! */
   this->layerStart = (LexNode **) New (&LexNet::_netHeap, this->nLayers * sizeof (LexNode *));

   this->nNodes = nn;
   this->node = (LexNode *) New (&LexNet::_netHeap, this->nNodes * sizeof (LexNode));

   /* pronId to pron mapping */ 
   this->pronlist = (Pron *) New (&LexNet::_netHeap, (this->voc->nprons + 1) * sizeof (Pron));
   this->pronlist[0] = NULL;

   /* initialise pointers to start of layers */
   LexNode *ln = this->node;
   for (int i = 0; i < this->nLayers; ++i) {
      layerCur[i] = ln;
      this->layerStart[i] = ln;
      ln += tnet->nNodesLayer[i];
   }

   /* initialise TLexNode -> LexNode pointers */
   ln = this->node;
   for (auto tln = tnet->root; tln; tln = tln->chain) {
      /* take node from appropriate layer memory block */
      tln->ln = layerCur[tln->layerId];
      ++layerCur[tln->layerId];
   }
   
   /* check whether all layer memory blocks are full (as they should be) */
   for (int i = 0; i < this->nLayers-1; ++i)
      assert (this->layerStart[i+1] == layerCur[i]);

   InitLMlaTree(this, tnet);

   /* the real conversion follows here: */
   this->start = tnet->start->ln;
   this->end = tnet->end->ln;
   
   int nlTotal = 0;
   /* copy info and convert links */
   for (auto tln = tnet->root; tln; tln = tln->chain) {
      ln = tln->ln;
      ln->type = (unsigned short) tln->type;

      switch (tln->type) {
	case LN_MODEL:
	  ln->data.hmm = tln->hmm;
	  break;
	case LN_CON:
	  ln->data.hmm = NULL;
	  break;
	case LN_WORDEND:
	  assert (tln->loWE == tln->hiWE);
	  ln->data.pron = tln->loWE;
	  this->pronlist[tln->loWE] = tln->pron;
	  break;
	default:
	  HError (9999, "HLVNet: unknown node type %d\n", tln->type);
	  break;
      }

      if (tnet->silDict) {
         this->lnSEsp = tnet->lnSEsp->ln;
         this->lnSEsil = tnet->lnSEsil->ln;
      }

      ln->lmlaIdx = tln->lmlaIdx;
      assert(ln->lmlaIdx < this->laTree->nNodes + this->laTree->nCompNodes);
      /*       printf ("LMLA idx %d lo %d hi %d\n", tln->lmlaIdx, tln->loWE, tln->hiWE); */

      /* only handle simple case here -- complex one has already been done in CreateCompLMLA() */
      if (ln->lmlaIdx < this->laTree->nNodes) {  
         this->laTree->node[ln->lmlaIdx].loWE = tln->loWE;
         this->laTree->node[ln->lmlaIdx].hiWE = tln->hiWE;
      }

      ln->nfoll = tln->nlinks;
      ln->foll = (LexNode **) New (&LexNet::_netHeap, ln->nfoll * sizeof (LexNode *));

      /* convert linked list of TLexLinks into array of pointers to
         end nodes of these links */
      int nl=0;
      LexNode **foll = ln->foll;
      for (auto tll = tln->link; tll; tll = tll->next, ++foll) {
         assert (tll->start == tln);
         *foll = tll->end->ln;
         ++nl;
      }
      assert (nl == ln->nfoll);
      nlTotal += nl;

      ++ln;
   }

   if (trace & T_NET)
      printf ("converted %d links\n", nlTotal);
}


TLexNet::TLexNet(Vocab *voc, HMMSet *hset, char *startWord, char *endWord, bool silDict):
  voc(voc), hset(hset), silDict(silDict) {

  startId = GetLabId (startWord, false);
  if (!startId) 
    HError (9999, "HLVNet: cannot find STARTWORD '%s'\n", startWord);

  endId = GetLabId (endWord, false);
  if (!endId) 
    HError (9999, "HLVNet: cannot find ENDWORD '%s'\n", endWord);
}

void TLexNet::init() {

  this->__init__();

  /* init phone sets A, B, AB, YZ and create LexNodes for AB and YZ */
  this->CollectPhoneStats ();

  /* Create initial phone (A) layer of z-a+b nodes also creates SA nodes*/
  this->CreateAnodes ();

  /* Create final phone (Z) layer of y-z+a nodes also creates ZS nodes */
  this->CreateZnodes ();

  /* Create silence (sil/sp) nodes connecting ZS and SA nodes */
  this->CreateSILnodes ();

  /* Create prefix tree nodes (B -- Y) */
  this->CreateBYnodes ();

  /* Create sentence start and end nodes */
  this->CreateStartEnd ();

  this->AssignWEIds();
}

void TLexNet::__init__() {

  heap = &TLexNet::tnetHeap;
  root = NULL;

  nlexA = nlexZ = nlexP = 0;
  nlexAB = nlexYZ = nlexZS = nlexSA = 0;
  nNodeA = nNodeZ = 0;

  lexA = (LabId *) New (&gcheap, LIST_BLOCKSIZE * sizeof(LabId));
  lexZ = (LabId *) New (&gcheap, LIST_BLOCKSIZE * sizeof(LabId));
  lexP = (LabId *) New (&gcheap, LIST_BLOCKSIZE * sizeof(LabId));

  for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
    lexABhash[i] = lexYZhash[i] = 
      lexZShash[i] = lexSAhash[i] = NULL;

  for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
    nodeAhash[i] = nodeZhash[i] = NULL;

  nNodeBY = 0;
  nodeBY = NULL;
  nNodeSIL = 0;
  nodeSIL = NULL;
  lmlaCount = 0;
  nPronIds = 0;

  for (int i = 0; i < NLAYERS; ++i)
    nNodesLayer[i] = 0;
}

void* TLexNet::operator new (size_t n) {
  return New (&TLexNet::tnetHeap, sizeof (TLexNet));
}

MemHeap TLexNet::tnetHeap("Net temp heap", MSTAK, 1, 0,100000, 800000);	/* used for temporary data in net creation */

void TLexNet::WriteTLex (char *fn)
{
   FILE *dotFile;
   bool isPipe;

   dotFile = FOpen (fn, NoOFilter, &isPipe);

   fprintf (dotFile, "digraph G {\n");
   fprintf (dotFile, "rankdir=LR;\n");
   fprintf (dotFile, "size=\"7,11\";\n");
   fprintf (dotFile, "{ A -> AB -> YZ -> Z -> ZS -> SA -> A2; }\n");

   /* A model nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "A;\n");
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
      for (auto ln = this->nodeAhash[i]; ln; ln = ln->next) {
         fprintf (dotFile, "n%p [label=\"%s\"];\n", ln, ln->lc ? ln->lc->name : "A?");
      }
   fprintf (dotFile, "}\n");

   /* links from A */
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
      for (auto ln = this->nodeAhash[i]; ln; ln = ln->next) {
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
         }
      }

   /* AB context nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "AB [shape=box];\n");
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexABhash[i]; ln; ln = ln->next)
         fprintf (dotFile, "n%p [label=\"%s-%s\", shape=box];\n", ln, ln->lc ? ln->lc->name : "A?", 
                  ln->rc ? ln->rc->name : "B?");
   fprintf (dotFile, "}\n");
   
   /* links from AB nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexABhash[i]; ln; ln = ln->next) {
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
         }
      }

   /* BY model nodes & links */
   for (auto ln = this->nodeBY; ln; ln = ln->next) {
      if (ln->type == LN_WORDEND) 
         fprintf (dotFile, "n%p [label=\"%s\", shape=diamond];\n", ln, ln->lc ? ln->lc->name : "BY?");
      else
         fprintf (dotFile, "n%p [label=\"%s\"];\n", ln, ln->lc ? ln->lc->name : "BY?");
      for (auto ll = ln->link; ll; ll = ll->next) {
         assert (ll->start == ln);
         fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
      }
   }


   /* YZ context nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "YZ [shape=box];\n");
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexYZhash[i]; ln; ln = ln->next)
         fprintf (dotFile, "n%p [label=\"%s-%s\", shape=box];\n", ln, ln->lc ? ln->lc->name : "Y?",
                  ln->rc ? ln->rc->name : "Z?");
   fprintf (dotFile, "}\n");


   /* links from YZ nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexYZhash[i]; ln; ln = ln->next) {
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
         }
      }

   /* Z model nodes & links */
   fprintf (dotFile, "Z;\n");
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
      for (auto ln = this->nodeZhash[i]; ln; ln = ln->next) {
         fprintf (dotFile, "n%p [label=\"%s\"];\n", ln, ln->lc ? ln->lc->name : "Z?");
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
         }
      }

   /* ZS context nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "ZS [shape=box];\n");
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexZShash[i]; ln; ln = ln->next)
         fprintf (dotFile, "n%p [label=\"%s-%s\", shape=box];\n", ln, ln->lc ? ln->lc->name : "Z?",
                  ln->rc ? ln->rc->name : "S?");
   fprintf (dotFile, "}\n");
   
   /* links from ZS nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexZShash[i]; ln; ln = ln->next)
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
         }

   /* SIL model nodes */
   fprintf (dotFile, "{ rank=same;\n");
   /*    fprintf (dotFile, "SIL ;\n"); */
   for (auto ln = this->nodeSIL; ln; ln = ln->next)
      fprintf (dotFile, "n%p [label=\"%s\"];\n", ln, ln->lc ? ln->lc->name : "SIL?");
   fprintf (dotFile, "}\n");

   /* links from SIL model nodes */
   for (auto ln = this->nodeSIL; ln; ln = ln->next)
      for (auto ll = ln->link; ll; ll = ll->next) {
         assert (ll->start == ln);
         fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
      }
   

   /* SA context nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "SA [shape=box];\n");
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexSAhash[i]; ln; ln = ln->next)
         fprintf (dotFile, "n%p [label=\"%s-%s\", shape=box];\n", ln, ln->lc ? ln->lc->name : "S?",
                  ln->rc ? ln->rc->name : "A?");
   fprintf (dotFile, "}\n");
   
   /* links from SA nodes */
   for (int i = 0; i < LEX_CON_HASH_SIZE; ++i)
      for (auto ln = this->lexSAhash[i]; ln; ln = ln->next)
         for (auto ll = ln->link; ll; ll = ll->next) {
            assert (ll->start == ln);
            if (ll->end->type == LN_WORDEND)
               fprintf (dotFile, "n%p -> n%p;\n", ll->start, ll->end);
            else
               fprintf (dotFile, "n%p -> n2%p;\n", ll->start, ll->end);
         }
   /* copy of A nodes */
   fprintf (dotFile, "{ rank=same;\n");
   fprintf (dotFile, "A2 [label=A];\n");
   for (int i = 0; i < LEX_MOD_HASH_SIZE; ++i)
      for (auto ln = this->nodeAhash[i]; ln; ln = ln->next) {
         fprintf (dotFile, "n2%p [label=\"%s\"];\n", ln, ln->lc ? ln->lc->name : "A?");
      }
   fprintf (dotFile, "}\n");

   fprintf (dotFile, "}\n");
   FClose (dotFile, FALSE);
}

/* -------------- vocab handling --------------- */

void MarkAllWordsfromLat (Vocab *voc, Lattice *lat, bool silDict)
{
   for (int i = 0; i < lat->nn; ++i) {
      LNode* ln = &lat->lnodes[i];
      Word word = ln->word;

      if (!word)
	continue;

      word->aux  = (Ptr) 1;
      if (silDict) { /* -/sp/sil style dict: only mark - variants */

	/* exclude !SENT_START/END */
	if ((word->nprons % 3) != 0)
	  continue;

	for (auto pron = word->pron; pron ; pron = pron->next) {
	  pron->aux  = (Ptr) 1; pron = pron->next;       /* sp variant */
	  pron->aux  = (Ptr) 0; pron = pron->next;       /* sil variant */
	  pron->aux  = (Ptr) 0;
	}
      }
      else {         /* lvx-style dict: mark all prons */
	for (auto pron = word->pron; pron ; pron = pron->next)
	  pron->aux  = (Ptr) 1;
      }
   }
}
