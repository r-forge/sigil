#!/usr/bin/env python
from extractor import Extractor
from input import InputReader, CQPFormat
import utils
import re

__VERSION__ = "$Id: biberextractor.py 149 2007-09-08 12:03:01Z sm $"

def compileRE(str):
        return re.compile(str, re.VERBOSE | re.IGNORECASE)
    

class BiberExtractor(Extractor):
    """Feature extractor for the features described in "Variation across
    speech and writing" by D. Biber. 
    
    The features are extracted from the words and POS tags using regular
    expressions.

    Usage: biberextractor.py [FILES]
    
    The [FILES] have to contain one token per line, with columns separated
    by a tab-character. The first column has to contain the word, the 
    second column the POS-tag (from the BNC C5 tagset), and the third
    column the lemmatized form.

    The output is written to a file named BiberExtractor_[UNIQUEID].tab,
    where [UNIQUEID] is determined from the file names given. This allows
    several extraction proccesses to be started on subsets of the data
    without interference.
    
    """
    # GENERAL PATTERNS
    #######################################################
    WORDPOS = r"((?<=\s)(%s)/(%s)(?=\s))"
    POS = WORDPOS % ("\S+","%s")
    WORD = WORDPOS % ("%s","\S{3}")
    
    # two words
    WORDS_2 = WORD + "\s" + WORD
    
    #three words
    WORDS_3 = WORDS_2 + "\s" + WORD
    
    WORDS_4 = WORDS_3 + "\s" + WORD
    
    # HELPER PATTERNS, e.g. simplified tags such as N 
    #######################################################
    # forms of "have", including ones with negative particle attached
    H_HAVE = (POS % "VH\w") + "(\s" + (POS % "XX0") + ")?" 
    
    # auxiliary verbs
    H_AUX = (POS % "V[^V].")
    
    # past tense or past participled verb
    H_VBN = (POS % "V[BDHV][DN]")
    
    # any form of "be"
    H_BE = (POS % "VB.")
    
    # any kind of punctuation or start/end of sentence
    H_PUN = "(</?s>|" + (POS % "PU.") + ")"
    
    # any word (but not a sentence boundary!)
    H_ANY = POS % "\S{3}" 
    
    
    # LISTS OF VERBS, taken from Quirk et al, 1985
    ######################################################
    # "public" verbs
    VERBS_PUBLIC  = frozenset( 
        ["acknowledge", "add", "admit", "affirm", "agree",
        "allege", "announce", "argue", "assert", "bet", "boast", "certify",
        "claim", "comment", "complain", "concede", "confess", "confide",
        "contend", "convey", "declare", "deny", "disclose", "exclaim",
        "explain", "forecast", "foretell", "guarantee", "hint", "insist",
        "maintain", "mention", "object", "predict", "proclaim", "promise",
        "pronounce", "prophesy", "protest", "remark", "repeat", "reply",
        "report", "retort", "say", "state", "submit", "suggest", "swear",
        "testify", "vow", "warn", "write"])
    
    # "private" verbs
    VERBS_PRIVATE = frozenset( 
        ["accept", "anticipate", "ascertain", "assume", "believe",
        "calculate", "check", "conclude", "conjecture", "consider", "decide",
        "deduce", "deem", "demonstrate", "determine", "discern", "discover",
        "doubt", "dream", "ensure", "establish", "estimate", "expect", "fancy",
        "fear", "feel", "find", "foresee", "forget", "gather", "guess", "hear",
        "hold", "hope", "imagine", "imply", "indicate", "infer", "insure",
        "judge", "know", "learn", "mean", "note", "notice", "observe",
        "perceive", "presume", "presuppose", "pretend", "prove", "realize",
        "reason", "recall", "reckon", "recognize", "reflect", "remember",
        "reveal", "see", "sense", "show", "signify", "suppose", "suspect",
        "think", "understand"])
    
    # suasive verbs
    VERBS_SUASIVE = frozenset(
        ["agree", "allow", "arrange", "ask", "beg", "command",
        "concede", "decide", "decree", "demand", "desire", "determine",
        "enjoin", "ensure", "grant", "insist", "instruct", "intend", "move",
        "ordain", "order", "pledge", "pray", "prefer", "pronounce", "propose",
        "recommend", "request", "require", "resolve", "rule", "stipulate",
        "suggest", "urge", "vote"])
    
    VERBS_SEEM = frozenset(["seem", "appear"])
    
    # Features as described in Appendix II of Biber, 1988.
    #######################################################
    
    #### TENSE AND ASPECT MARKRES ####
    
    # 01  past tense, identified by POS tags
    RE_PAST_TENSE = POS % "V[BDHV]D" 
    
    # 02 perfect aspect
    # HAVE + (ADV) + (ADV) + VBN
    # HAVE + N/PRO + VBN
    RE_PERFECT_ASPECT_1 = H_HAVE + \
                          "(\s" +  (POS % "AV\w") + "){0,2}" + \
                          "\s" + H_VBN
    RE_PERFECT_ASPECT_2 = H_HAVE + \
                          "\s" + \
                          (POS % "NN\w|PN\w") + \
                          "\s" + \
                          H_VBN
    RE_PERFECT_ASPECT = "(" + RE_PERFECT_ASPECT_1 + "|" + \
                              RE_PERFECT_ASPECT_2 + ")"
    
    # 03 present tense
    RE_PRESENT_TENSE = POS % "V[BDHV][BZ]" 
    
    
    #### PLACE AND TIME ADVERBIALS ####
    
    # 04 place adverbials (not restricted to adverb POS)
    RE_PLACE_ADVERBIALS = WORD % '|'.join([
        'aboard', 'above', 'abroad', 'across', 'ahead', 'alongside', 'around',
        'ashore', 'astern', 'away', 'behind', 'below', 'beneath', 'beside',
        'downhill', 'downstairs', 'downstream', 'east', 'far', 'hereabouts', 
        'indoors', 'inland', 'inshore', 'inside', 'locally', 'near', 'nearby',
        'north', 'nowhere', 'outdoors', 'outside', 'overboard',  'overland',
        'overseas' 'south' 'underfoot', 'underground', 'underneath', 'uphill',
        'upstairs' 'upstream' 'west'])
    
    # 05 time adverbials (not restricted to adverb POS)
    RE_TIME_ADVERBIALS = WORD % '|'.join([
        'afterwards', 'again', 'earlier', 'early', 'eventually', 'formerly', 
        'immediately', 'initially', 'instantly', 'late', 'lately', 'later',
        'momentarily', 'now', 'nowadays', 'once', 'originally', 'persently',
        'previously', 'recently', 'shortly', 'simultaneously', 'soon', 
        'subsequently', 'today', 'tomorrow', 'tonight', 'yesterday'])
    
    
    #### PRONOUNS AND PRO-VERBS ####
    
    # 06 first person pronouns 
    RE_FIRST_PERSON_PRONOUNS = WORDPOS % ('|'.join([
        'I', 'me', 'we', 'us', 'my', 'our', 'myself', 'ourselves']),
        "PNP|PNX|DPS")
    
    # 07 second person pronouns 
    RE_SECOND_PERSON_PRONOUNS = WORDPOS % ('|'.join([
        'you', 'your', 'yourself', 'yourselves']), 
        "PNP|PNX|DPS")
    
    # 08 third person personal pronouns 
    RE_THIRD_PERSON_PRONOUNS = WORDPOS % ('|'.join([
        'she', 'he', 'they', 'her', 'him' 'them', 'his', 'their', 'himself',
        'herself', 'themselves']), 
        "PNP|PNX|DPS")
    
    # 09 pronoun "it"
    RE_PRONOUN_IT = WORDPOS % ('it', "PNP")
    
    # 10 demonstrative pronouns
    # POS restricted to DT0 to exclude "that" as CJT
    RE_DEMONSTRATIVE_PRONOUN = WORDPOS % ("that|this|these|those","DT0") + \
                                 "\s" + (POS % "V..|PU[NR]|PNQ|CJC")
    
    # 11 indefinite pronouns
    RE_INDEFINITE_PRONOUN = POS % "PNI"
    # Could also have used the following list given by Biber
    # ['anybody', 'anyone', 'anything', 'everybody', 'everyone', 'everything',
    # 'nobody', 'none', 'nothing', 'nowhere', 'somebody', 'someone', 
    # 'something'])
    
    # 12 pro-verb do
    # "do" when not as auxialiary or question
    RE_PROVERB_DO = "(?<!" + \
                        "(PU.|PNQ)" + \
                    "\s)" + \
                    (POS % "VD.") + \
                    "(?!\s" + \
                        "(" + (POS % "AV.") + "\s)?" + \
                        (POS % "V..")+ \
                    ")"
    
    #### QUESTIONS #### 
    
    # 13 direct WH-questions
    RE_WH_QUESTION = "<s>\s" + (POS % "PNQ|AVQ") + "\s" + H_AUX
    
    
    #### NOMINAL FORMS ####
    
    # 14 nominalizations
    RE_NOMINALIZATION = WORDPOS % ("\w+(" + '|'.join([
        "tions?","ments?","ness(es)?","it(y|ies)"]) + ")", "N..")
    
    # 15 gerunds
    RE_GERUNDS = "(?<!VB.\s)" + POS % ("VVG")
    
    
    # 16 all other nouns 
    # count all nouns here and then subtract later 
    RE_NOUNS = POS % "N.."
    
    
    #### PASSIVES ####
    
    RE_PASSIVES_1 = H_BE + "\s" + \
                              "(" + (POS % "AV.") + "\s)*" + \
                              H_VBN
    
    RE_PASSIVES_2 = H_BE + "\s"\
                              "(" + (POS % "N..") + "|" + (POS % "PN.") +")" +\
                              "\s" + H_VBN
                              
    RE_PASSIVES = "(" + RE_PASSIVES_1 + "|" + RE_PASSIVES_2 + ")"
    
    # 17 agentless passives
    RE_AGENTLESS_PASSIVES = RE_PASSIVES + "(?!\s" + (WORD % "by") + ")"
    
    # 18 "by" passives
    RE_BY_PASSIVES = RE_PASSIVES + "\s" + (WORD % "by")
    
    
    #### STATIVE FORMS ####
    
    # 19 "be" as main verb
    # this misses some forms ...  
    RE_BE_MAIN_VERB = (POS % "VB.") + "\s" + \
            (POS % "AT0|DPS|PR[FP]|AJ[0CS]|NP0")
                                 
    # 20 existential there
    RE_EXISTENTIAL_THERE = POS % "EX0"
    
    
    #### SUBORDINATION #### 
    
    # 21 "that" verb complements
    RE_THAT_VERB_COMP_1 = (POS % "PU.|CJC") + "\s" + \
                          (WORD % "that") + "\s" + \
                          (POS % "AT0|DT.|PN[IP]|EX0|NN[02]|NP0")
    
    RE_THAT_VERB_COMP_2 = (POS % "VV.") + "\s" + \
                          (WORD % "that") + "\s" + \
                          "(?!" + \
                              (POS % "V..|PUN|CJC") + \
                          ")"
                          
    RE_THAT_VERB_COMP_3 = (POS % "VV.") + "\s" + \
                          (POS % "PR[PF]") + "\s" + \
                          "(" + (POS % "[^N]..") + "\s)*?" + \
                          (WORDPOS % ("that","CJT"))
                          
    
    RE_THAT_VERB_COMP = "(" + RE_THAT_VERB_COMP_1 + "|" + \
                              RE_THAT_VERB_COMP_2 + "|" + \
                              RE_THAT_VERB_COMP_3 + ")"
    
    # 22 "that" adjective complements
    RE_THAT_ADJ_COMP = (POS % "AJ.") + "\s" + \
                       (WORD % "that")
    
    # 23 WH-clauses
    # Wh-determiners (DTQ) are left out
    RE_WH_CLAUSE = (POS % "VV.") + "\s" + \
                   (POS % "PNQ|AVQ") + "\s" + \
                   "(?!" + H_AUX + ")"
    
    # 24 infinitives
    RE_INFINITIVES = (POS % "TO0") + "\s" + \
                     "(" + (POS % "AV.") + "\s)*" + (POS % "V.I")
    
    # 25 present participial clauses
    RE_PRESENT_PARTICIPLE = H_PUN + "\s" + \
                            (POS % "V.G") + "\s" + \
                            (POS % "PR.|AT0|DT.|DTQ|PN.|AV.")
    
    # 26 past participial clauses
    RE_PAST_PARTICIPLE = H_PUN + "\s" + \
                         (POS % "V.N") + "\s" + \
                         (POS % "PR.|AV.")
    
    # 27 past participial WHIZ deletion relatives
    RE_PAST_PARTICIPLE_WHIZ = (POS % "N..|PNI") + "\s" + \
                              (POS % "V.N") + "\s" + \
                              (POS % "PR.|VB.|AV.")
    
    # 28 present participial WHIZ deletion relatives
    RE_PRESENT_PARTICIPLE_WHIZ = (POS % "N..|PNI") + "\s" + \
                              (POS % "V.G")
                              
    # 29 "that" relative clauses on subject position
    RE_THAT_SUBJ =  (POS % "N..") + "\s" + \
                    (WORD % "that") + "\s" + \
                    "("+ (POS % "AV.") +  "\s)*" + \
                    (POS % "V..")
    
    # 30 "that" relative clauses on object position
    RE_THAT_OBJ =   (POS % "N..") + "\s" + \
                    (WORD % "that") + "\s" + \
                    (POS % "AT0|DT.|PNP|NN[02]|NP0|AJ.")
    
    # 31 WH relative clauses on subject position
    # ignoring some of the restrictions made by Biber
    RE_WH_SUBJ = H_ANY + "\s" + H_ANY + "\s" + \
                 (POS % "N..") + "\s" + \
                 (POS % "PNQ|DTQ") + "\s" + \
                 "("+ (POS % "AV.") +  "\s)*" + \
                 (POS % "V..")
    
    # 32 WH relative clauses on object position
    # ignoring some of the restrictions made by Biber
    RE_WH_OBJ  = H_ANY + "\s" + H_ANY + "\s" + \
                 (POS % "N..") + "\s" + \
                 (POS % "PNQ|DTQ") + "\s" + \
                 "(?!" + (POS % "V..|AV.")+ ")"
    
    # 33 pied-piping relatives 
    RE_PIED_PIPING = (POS % "PR.") + "\s" + (POS % "PNQ|DTQ")
    
    # 34 sentence relatives 
    RE_SENTENCE_RELATIVES = (WORD % ",") + "\s" + (POS % "PNQ|DTQ")
    
    
    #### ADVERBIAL CLAUSES ####
    
    # 35 causative adverbial subordinator: "because" 
    RE_BECAUSE = WORD % "because"
    
    # 36 concessive adverbial subordinators: "although", "though"
    RE_THOUGH = WORD % "(al)?though"
    
    # 37 conditional adverbial subordinators
    RE_IF = WORD % "if|unless"
    
    # 38 other adverbial subordinators
    RE_OTHER_ADV_SUB = "(" + \
        (WORDPOS % ("since|while|whilst|whereas|whereby","CJS")) + "|" + \
        "(" + '|'.join([
                    WORDS_2 % ("such","that"),
                    WORDS_2 % ("so","that"),
                    WORDS_2 % ("inasmuch","as"),
                    WORDS_2 % ("forasmuch","as"),
                    WORDS_2 % ("insofar","as"),
                    WORDS_2 % ("insomuch","as"),
                    WORDS_3 % ("as","long","as"),
                    WORDS_3 % ("as","soon","as")
                ]) + ")" + ")"
    
    
    #### PREPOSITIONAL PHRASES #### 
    
    # 39 total prepositional phrases
    RE_PREP = POS % "PR."
    
    
    #### ADJECTIVES AND ADVERBS #### 
    
    # 40 attributive adjectives 
    # all non-predicative adjectives, calculated later
    
    # 41 predicative adjectives 
    RE_ADJ_PRED = H_BE + "\s" + \
                  (POS % "AJ.") + "\s" + \
                  "(" + (POS % "AV.")+"\s)*" + \
                  "(?!" + (POS % "AJ.|N..")+")"
    
    # total adjectives 
    RE_ADJ = POS % "AJ."
    
    # 42 total adverbs
    RE_ADV = POS % "AV."
    
    #### LEXICAL SPECIFICITY #### 
    # 43 type/token ratio
    # ... calculated elsewhere
    # 44 word length
    # ... also
    
    #### LEXICAL CLASSES #### 
    
    # 45 conjuncts 
    
    RE_CONJUNCTS_1 = WORD % '|'.join([
        'alternatively', 'altogether', 'consequently', 'conversely', 'eg',
        r'e\.g\.', 'else', 'furthermore', 'hence', 'however', r'i\.e\.',
        'instead', 'likewise', 'moreover', 'namely', 'nevertheless', 
        'nonetheless', 'notwithstanding', 'otherwise', 'rather', 'similarly',
        'therefore', 'thus', 'viz'])
    
    RE_CONJUNCTS_2 = "(" + '|'.join([
        WORDS_2 % ("in","comparison|contrast|particular|addition" + \
                   "|conclusion|consequence|sum|summary"), 
        WORDS_3 % ("in","any","case|event"),
        WORDS_3 % ("in", "other", "words"),
        WORDS_2 % ("for", "example|instance"), 
        WORDS_2 % ("by", "contrast|comparison"),
        WORDS_3 % ("as", "a", "result|consequence"),
        WORDS_3 % ("on", "the", "contrary"),
        WORDS_4 % ("on","the","other","hand") ]) + ")"
        # some left out here ... 
    
    RE_CONJUNCTS = H_PUN + "\s" + \
                   "(" + RE_CONJUNCTS_1 + "|" + RE_CONJUNCTS_2 + ")"
    
    # 46 downtoners 
    RE_DOWNTONERS = WORDPOS % ("almost|barely|hardly|merely|mildly|nearly" + \
                           "|only|partially|partly|practically|scarcely" + \
                           "|slightly|somewhat", "AV0")
                           
    # 47 hedges
    RE_HEDGES = "(" + '|'.join([
        WORDS_2 % ("at","about"),
        WORDS_2 % ("something","like"),
        WORDS_3 % ("more","or","less"),
        WORD % "almost|maybe",
        "(" + "(?<!(AT0|DT0|AJ.|PNP|AVQ)\s)" + WORDS_2 % ("sort|kind", "of") + ")" 
        ]) + ")"
    
    # 48 amplifiers 
    RE_AMPLIFIERS = WORD % ("absolutely|altogether|completely|enormously" + \
                           "|entirely|extremely|fully|greatly|highly" + \
                           "|intensely|perfectly|strongly|thoroughly" + \
                           "|totally|utterly|very")
    
    # 49 emphatics 
    RE_EMPHATICS = "(" + '|'.join([
        WORDS_2 % ("for", "sure"),
        WORDS_2 % ("a", "lot"),
        WORD % ("just|really|most|more"),
        "(" + (POS % "VD.") + "\s" + (POS % ("VV.")) +")",
        "(" + (WORD % "real|so") + "\s" + (POS % "AJ.")+ ")" ]) + ")"
    
    # 50 discourse particles 
    RE_DISCOURSE_PARTICLES = H_PUN + "\s" + (WORD % "well|now|anyway|anyhow|anyways")
    
    # 51 demonstratives (but not demonstrative pronouns)
    RE_DEMONSTRATIVES = WORDPOS % ("that|this|these|those","DT0") + \
                                 "(?!\s" + (POS % "V..|PU[NR]|PNQ|CJC") + ")"
    
    
    #### MODALS ####
    
    # 52 possibility modals 
    RE_MODAL_POSSIBILITY = WORDPOS % ("can|may|might|could", "VM0")
    
    # 53 necessity modals
    RE_MODAL_NECESSITY = WORDPOS % ("ought|should|must", "VM0")
    
    # 54 predictive modals 
    RE_MODAL_PREDICTIVE = WORDPOS % ("will|would|shall|'ll","VM0")
    
    
    #### SPECIALIZED VERB CLASSES #### 
    # see above for wordlists, extracted from lemmatized form elsewhere
    # 55 public verbs
    # 56 private verbs
    # 57 suasive verbs 
    # 58 seam/appear
    
    
    #### REDUCED FORMS #### 
    
    # 59 contractions
    RE_CONTRACTIONS = WORDPOS % ("\w*'\w+","PN.|XX0|V[^V].")
    
    # 60 subordinator-"that" deletion
    # left out
    
    # 61 stranded preposistions
    RE_STRANDED_PREPOSITION = POS % ("PR.") + "\s" + "(</s>|" + (POS % "PU[NR]") + ")"
    
    # 62 split infinitives
    RE_SPLIT_INFINITIVE = POS % ("TO0") + "\s" + \
                          "(" + POS % ("AV.") + "\s)+" + \
                          POS % ("VV[IB]")  
    
    # 63 split auxiliary
    RE_SPLIT_AUXILIARY = POS % ("V[^V].") + "\s" + \
                          "(" + POS % ("AV.") + "\s)+" + \
                          POS % ("VV[IB]")  
    
    # 64 phrasal coordination
    RE_PHRASAL_COORDINATION = r"\s \S+/(AV|NN|NP|AJ|V)..? \s \S+/CJC \s \S+/\1..? \s"
    
    # 65 independent clause coordination
    # left out
    
    # 66 synthetic negation
    RE_NEG_SYNTHETIC = "(" + (WORD % "no") + "\s" + (POS % "AJ.|N..") + "|" + \
                             (WORD % "neither|nor") + ")"
    
    # 66 analytic negation 
    RE_NEG_ANALYTIC = POS % "XX0"
    
    # List of all features
    DIRECT_FEATS = [ 
             ("01_past_tense", RE_PAST_TENSE), 
             ("02_perfect_aspect", RE_PERFECT_ASPECT),
             ("03_present_tense", RE_PRESENT_TENSE),
             ("04_place_adverbials", RE_PLACE_ADVERBIALS),
             ("05_time_adverbials", RE_TIME_ADVERBIALS),
             ("06_first_person_pronouns", RE_FIRST_PERSON_PRONOUNS),
             ("07_second_person_pronouns", RE_SECOND_PERSON_PRONOUNS),
             ("08_third_person_pronouns", RE_THIRD_PERSON_PRONOUNS),
             ("09_pronoun_it", RE_PRONOUN_IT),
             ("10_demonstrative_pronoun", RE_DEMONSTRATIVE_PRONOUN),
             ("11_indefinite_pronoun", RE_INDEFINITE_PRONOUN),
             ("12_proverb_do", RE_PROVERB_DO),
             ("13_wh_question", RE_WH_QUESTION),
             ("14_nominalization", RE_NOMINALIZATION),
             ("15_gerunds", RE_GERUNDS),
             ("17_agentless_passives", RE_AGENTLESS_PASSIVES), 
             ("18_by_passives", RE_BY_PASSIVES),
             ("19_be_main_verb", RE_BE_MAIN_VERB),
             ("20_existential_there", RE_EXISTENTIAL_THERE),
             ("21_that_verb_comp", RE_THAT_VERB_COMP),
             ("22_that_adj_comp",RE_THAT_ADJ_COMP),
             ("23_wh_clause", RE_WH_CLAUSE),
             ("24_infinitives", RE_INFINITIVES),
             ("25_present_participle", RE_PRESENT_PARTICIPLE),
             ("26_past_participle", RE_PAST_PARTICIPLE),
             ("27_past_participle_whiz", RE_PAST_PARTICIPLE_WHIZ),
             ("28_present_participle_whiz", RE_PRESENT_PARTICIPLE_WHIZ),
             ("29_that_subj", RE_THAT_SUBJ), 
             ("30_that_obj", RE_THAT_OBJ),
             ("31_wh_subj", RE_WH_SUBJ),
             ("32_wh_obj", RE_WH_OBJ),
             ("33_pied_piping", RE_PIED_PIPING),
             ("34_sentence_relatives", RE_SENTENCE_RELATIVES),
             ("35_because", RE_BECAUSE),
             ("36_though", RE_THOUGH),
             ("37_if", RE_IF), 
             ("38_other_adv_sub", RE_OTHER_ADV_SUB),
             ("39_prepositions", RE_PREP), 
             ("41_adj_pred", RE_ADJ_PRED),
             ("42_adverbs", RE_ADV),
             ("45_conjuncts", RE_CONJUNCTS),
             ("46_downtoners", RE_DOWNTONERS),
             ("47_hedges", RE_HEDGES),
             ("48_amplifiers", RE_AMPLIFIERS),
             ("49_emphatics", RE_EMPHATICS),
             ("50_discourse_particles", RE_DISCOURSE_PARTICLES),
             ("51_demonstratives", RE_DEMONSTRATIVES), 
             ("52_modal_possibility", RE_MODAL_POSSIBILITY),
             ("53_modal_necessity", RE_MODAL_NECESSITY),
             ("54_modal_predictive", RE_MODAL_PREDICTIVE), 
             ("59_contractions", RE_CONTRACTIONS), 
             ("61_stranded_preposition", RE_STRANDED_PREPOSITION),
             ("62_split_infinitve", RE_SPLIT_INFINITIVE),
             ("63_split_auxiliary", RE_SPLIT_AUXILIARY), 
             ("64_phrasal_coordination", RE_PHRASAL_COORDINATION), 
             ("66_neg_synthetic", RE_NEG_SYNTHETIC), 
             ("67_neg_analytic", RE_NEG_ANALYTIC)
             ]
    
    ## Pre-compile REs to speed things up a little
    DIRECT_FEATS = [(name, compileRE(regex)) for (name,regex) in DIRECT_FEATS]
    
    # features that need to be counted for calculations but are not directly
    # used.
    CALC_FEATS = [
        ("nouns", RE_NOUNS),
        ("adjectives", RE_ADJ)
        ] 
   
    ## Pre-compile REs to speed things up a little
    CALC_FEATS = [(name, compileRE(regex)) for (name,regex) in CALC_FEATS]
   
             
    # the names of the features that are calculated on the basis of other
    # features 
    CALCULATED_FEATS = ["16_other_nouns", "40_adj_attr", "43_type_token", "44_mean_word_length"]
    
    # features to be counted on lemmatized from
    # list of tuples of the form (featName,wordset) 
    LEMMA_FEATS = [("55_verb_public", VERBS_PUBLIC), 
                   ("56_verb_private", VERBS_PRIVATE), 
                   ("57_verb_suasive", VERBS_SUASIVE), 
                   ("58_verb_seem", VERBS_SEEM)]
    
    
    ### further possible non-Biber features
    # - all past participles (V.N)
    # - total subordinators (CJS) 
    
        
    def getName(self):
        return "BiberExtractor"
    
    def getFeatureNames(self):
        return self.__featureNames
    
    def getAvailableOptions(self):
        short = ""
        long = []
        return short, long
        
    def __init__(self):
        self.__featureNames = sorted(
            [name for (name,re) in self.DIRECT_FEATS]
          + [name for (name,re) in self.LEMMA_FEATS]  
          + self.CALCULATED_FEATS
            )
        Extractor.__init__(self)
        
    def disambiguatePOS(self,poslist):
        # For ambigious tags, use the more probable analysis (the first tag)
        return [x.split('-')[0] for x in poslist]
    
    def getWordsWithPOS(self,wordlist,poslist):
        if len(wordlist) != len(poslist):
            raise Error("Word- and POS-list are not of the same length!")
        return [word + "/" + pos for (word,pos) in zip(wordlist,poslist)]
    
    def extractWithREs(self,features,text):
        featDict = {}
        for (featName, featRE) in features:
            #print featName
            #print featRE.pattern
            matches = featRE.findall(text)
            #matches = re.findall("(" + featRE + ")", text, re.VERBOSE | re.IGNORECASE)
            #self.showMatches(featRE,text)
            count = len(matches)
            featDict[featName] = count
        return featDict
    
    def extractFromLemmatatizedForms(self,features,lemmaList):
        featDict = {}
        for (featName, words) in features:
            count = sum([1 for lemma in lemmaList if lemma in words])
            featDict[featName] = count
        return featDict
    
    def extractStatistics(self,cqpf):
        featDict = {}
        wordList = cqpf.getColumn(0)
        posList = cqpf.getColumn(1)
        lemmaList = cqpf.getColumn(2) # use lemmatized forms as tokens
        tokens = len(lemmaList)
        types = len(frozenset(lemmaList))
        featDict["43_type_token"] = types / float(tokens)
        wordLengths = []
        for i in range(len(wordList)):
            if posList[i] not in ["PUN", "PUL", "PUR", "PUQ"]:
                wordLengths.append(len(wordList[i]))
        featDict["44_mean_word_length"] = utils.mean(wordLengths)
        
        return featDict
    
    
    
    # calculate some features after extraction
    # the passed dictionary will be changed!
    def calculateFeats(self, featDict):
        featDict["16_other_nouns"] = (
            featDict["nouns"] 
            - featDict["14_nominalization"] 
            - featDict["15_gerunds"])
        del featDict["nouns"]
        
        featDict["40_adj_attr"] = (
            featDict["adjectives"] 
            - featDict["41_adj_pred"]  )
        del featDict["adjectives"]
        
    def normalizeByLength(self,featDict,length):
        for feat in featDict:
            featDict[feat] = featDict[feat] / float(length) * 1000 
    
    def showMatches(self,pattern,text,num=10, context = 40):
        count = 0
        for match in pattern.finditer(text):
            start = match.start()
            end = match.end()
            print text[start-context:start] + " [ " + text[start:end] + " ] " + text[end:end+context]
            count += 1
            if count == num:
                return
    
    def process(self,file):
        feats = {}
        Extractor.process(self,file)
        ir = InputReader(file)
        ir.read()
        cqpf = CQPFormat(ir.getText())
        #words = ' '.join(cqpf.getColumn(0))
        #pos = ' '.join(self.disambiguatePOS(cqpf.getColumn(1)))
        lemma = cqpf.getColumn(2)
        sentences = cqpf.getAnnotations("s")
        wordpostmp = []
        for (start,end,attr) in sentences:
            wordpostmp.append('<s>')
            wordpostmp.extend(self.getWordsWithPOS(
                                cqpf.getColumn(0)[start:end],
                                self.disambiguatePOS(cqpf.getColumn(1)[start:end])))
            wordpostmp.append('</s> ')
        wordpos = ' '.join(wordpostmp)
        feats.update(self.extractWithREs(self.DIRECT_FEATS,wordpos))
        feats.update(self.extractWithREs(self.CALC_FEATS,wordpos))
        feats.update(self.extractFromLemmatatizedForms(self.LEMMA_FEATS,lemma))
        self.calculateFeats(feats)
        self.normalizeByLength(feats, len(lemma))
        feats.update(self.extractStatistics(cqpf))
        # print feats
        return ir.getID(),feats
        
if __name__ == "__main__":
    BiberExtractor()
