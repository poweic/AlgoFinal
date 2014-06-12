class HiddenMarkovModel {
  public:
    HiddenMarkovModel(char* hmmList_fn, char* mmf_fn, char* inXFormDir_fn,
	char* hmmDir, char* hmmExt) {

      /* init model heap & set early to support loading MMFs */
      CreateHeap(&_modelHeap, "Model heap",  MSTAK, 1, 0.0, 100000, 800000 );
      CreateHMMSet(&_hset, &_modelHeap, TRUE); 

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
      ::ConvDiagC (&_hset, TRUE);
      ::ConvLogWt (&_hset);

      if (trace&T_TOP)
	printf("Read %d physical / %d logical HMMs\n", _hset.numPhyHMM, _hset.numLogHMM);  
    }

    void SetStreamWidths(Boolean *eSep) {
      ::SetStreamWidths(_hset.pkind, _hset.vecSize, _hset.swidth, eSep);
    }

    Boolean UpdateSpkrStats(XFInfo *xfinfo, char *datafn) {
      return ::UpdateSpkrStats(&_hset, xfinfo, datafn);
    }

    MLink FindMacroStruct(char type, Ptr structure) {
      return ::FindMacroStruct (&_hset, type, structure);
    }

    MLink FindMacroName(char type, LabId id) {
      return ::FindMacroName(&_hset, type, id);
    }

    HMMSet& get_hset() { return _hset; }
  private:
    MemHeap _modelHeap;
    HMMSet _hset;
};

class Vocabulary {
  public:
    Vocabulary(char* dict_fn) {
      this->init(dict_fn);
    }

    void init(char* dict_fn) {

      /* Read dictionary */
      if (trace & T_TOP)
	printf ("Reading dictionary from %s\n", dict_fn);

      InitVocab (&_vocab);
      if (ReadDict (dict_fn, &_vocab) < SUCCESS)
	HError (9999, "Initialise: ReadDict failed");
    }

    void ConvertSilDict (LabId spLab, LabId silLab, LabId startLab, LabId endLab) {
      ::ConvertSilDict(&_vocab, spLab, silLab, startLab, endLab);
    }

    void MarkAllProns () {
      ::MarkAllProns(&_vocab);
    }

    void MarkAllWords () {
      ::MarkAllWords(&_vocab);
    }

    void UnMarkAllWords () {
      ::UnMarkAllWords(&_vocab);
    }

    void MarkAllWordsfromLat (Lattice* lat, Boolean silDict) {
      ::MarkAllWordsfromLat(&_vocab, lat, silDict);
    }

    Vocab& get_vocab() {
      return _vocab;
    }

  private:
    Vocab _vocab;
};

class LanguageModel {
  public:
    LanguageModel(Vocabulary& vocab): lm(NULL), _vocab(vocab) {
      CreateHeap (&_lmHeap, "LM heap", MSTAK, 1, 0,1000000, 10000000);
    }

    void loadFromFile(char* lm_fn) {
      if (!lm_fn) 
	HError (9999, "HDecode: no LM or lattice specified");
      lm = CreateLM (&_lmHeap, lm_fn, startWord, endWord, &_vocab.get_vocab());
    }

    void loadFromLattice(char* lat_fn, Lattice* lat) {
      lm = CreateLMfromLat (&_lmHeap, lat_fn, lat, &_vocab.get_vocab());
    }

    FSLM* get_lm() {
      return lm;
    }

    Lattice *ReadLattice(FILE *file, Boolean shortArc, Boolean add2Dict) {
      return ::ReadLattice(file, &_lmHeap, &_vocab.get_vocab(), shortArc, add2Dict);
    }

    void ResetHeap() {
      ::ResetHeap(&_lmHeap);
    }

  private:
    FSLM *lm;         /* language model */
    MemHeap _lmHeap;

    Vocabulary& _vocab;
};

class Decoder {
  public:
    Decoder(LanguageModel& lm, Vocabulary& vocab, HiddenMarkovModel& hmm);
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

    DecoderInst* _decInst;

    MemHeap _netHeap;
    // MemHeap _lmHeap;
    MemHeap _inputBufHeap;
    MemHeap _transHeap;
    MemHeap _regHeap;
};

