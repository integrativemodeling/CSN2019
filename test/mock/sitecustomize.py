import ihm
try:
    import ihm.reference
except ImportError:
    pass
import os

# If we're running from an SGE job, override the from_pubmed_id() and
# from_accession() functions to return cached values, since we don't have
# network access (needed to query PubMed or UniProt directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
        pmid=32034103,
        title='Structural dynamics of the human COP9 signalosome revealed '
              'by cross-linking mass spectrometry and integrative modeling.',
        journal='Proc Natl Acad Sci U S A', volume=117,
        page_range=('4088', '4098'), year=2020,
        authors=['Gutierrez C', 'Chemmama IE', 'Mao H', 'Yu C', 'Echeverria I',
                 'Block SA', 'Rychnovsky SD', 'Zheng N', 'Sali A', 'Huang L'],
        doi='10.1073/pnas.1915542117')

# Auto-generated by util/get_accession_map.py
accession_map = {
    'Q13098': {'db_code':'CSN1_HUMAN', 'sequence':
        'MPLPVQVFNLQGAVEPMQIDVDPQEDPQNAPDVNYVVENPSLDLEQYAASYSGLMRIERLQFIADH'
        'CPTLRVEALKMALSFVQRTFNVDMYEEIHRKLSEATRSSLRELQNAPDAIPESGVEPPALDTAWVE'
        'ATRKKALLKLEKLDTDLKNYKGNSIKESIRRGHDDLGDHYLDCGDLSNALKCYSRARDYCTSAKHV'
        'INMCLNVIKVSVYLQNWSHVLSYVSKAESTPEIAEQRGERDSQTQAILTKLKCAAGLAELAARKYK'
        'QAAKCLLLASFDHCDFPELLSPSNVAIYGGLCALATFDRQELQRNVISSSSFKLFLELEPQVRDII'
        'FKFYESKYASCLKMLDEMKDNLLLDMYLAPHVRTLYTQIRNRALIQYFSPYVSADMHRMAAAFNTT'
        'VAALEDELTQLILEGLISARVDSHSKILYARDVDQRSTTFEKSLLMGKEFQRRAKAMMLRAAVLRN'
        'QIHVKSPPREGSQGELTPANSQSRMSTNM'},
    'P61201': {'db_code':'CSN2_HUMAN', 'sequence':
        'MSDMEDDFMCDDEEDYDLEYSEDSNSEPNVDLENQYYNSKALKEDDPKAALSSFQKVLELEGEKGE'
        'WGFKALKQMIKINFKLTNFPEMMNRYKQLLTYIRSAVTRNYSEKSINSILDYISTSKQMDLLQEFY'
        'ETTLEALKDAKNDRLWFKTNTKLGKLYLEREEYGKLQKILRQLHQSCQTDDGEDDLKKGTQLLEIY'
        'ALEIQMYTAQKNNKKLKALYEQSLHIKSAIPHPLIMGVIRECGGKMHLREGEFEKAHTDFFEAFKN'
        'YDESGSPRRTTCLKYLVLANMLMKSGINPFDSQEAKPYKNDPEILAMTNLVSAYQNNDITEFEKIL'
        'KTNHSNIMDDPFIREHIEELLRNIRTQVLIKLIKPYTRIHIPFISKELNIDVADVESLLVQCILDN'
        'TIHGRIDQVNQLLELDHQKRGGARYTALDKWTNQLNSLNQAVVSKLA'},
    'Q9UNS2': {'db_code':'CSN3_HUMAN', 'sequence':
        'MASALEQFVNSVRQLSAQGQMTQLCELINKSGELLAKNLSHLDTVLGALDVQEHSLGVLAVLFVKF'
        'SMPSVPDFETLFSQVQLFISTCNGEHIRYATDTFAGLCHQLTNALVERKQPLRGIGILKQAIDKMQ'
        'MNTNQLTSIHADLCQLCLLAKCFKPALPYLDVDMMDICKENGAYDAKHFLCYYYYGGMIYTGLKNF'
        'ERALYFYEQAITTPAMAVSHIMLESYKKYILVSLILLGKVQQLPKYTSQIVGRFIKPLSNAYHELA'
        'QVYSTNNPSELRNLVNKHSETFTRDNNMGLVKQCLSSLYKKNIQRLTKTFLTLSLQDMASRVQLSG'
        'PQEAEKYVLHMIEDGEIFASINQKDGMVSFHDNPEKYNNPAMLHNIDQEMLKCIELDERLKAMDQE'
        'ITVNPQFVQKSMGSQEDDSGNKPSSYS'},
    'Q9BT78': {'db_code':'CSN4_HUMAN', 'sequence':
        'MAAAVRQDLAQLMNSSGSHKDLAGKYRQILEKAIQLSGAEQLEALKAFVEAMVNENVSLVISRQLL'
        'TDFCTHLPNLPDSTAKEIYHFTLEKIQPRVISFEEQVASIRQHLASIYEKEEDWRNAAQVLVGIPL'
        'ETGQKQYNVDYKLETYLKIARLYLEDDDPVQAEAYINRASLLQNESTNEQLQIHYKVCYARVLDYR'
        'RKFIEAAQRYNELSYKTIVHESERLEALKHALHCTILASAGQQRSRMLATLFKDERCQQLAAYGIL'
        'EKMYLDRIIRGNQLQEFAAMLMPHQKATTADGSSILDRAVIEHNLLSASKLYNNITFEELGALLEI'
        'PAAKAEKIASQMITEGRMNGFIDQIDGIVHFETREALPTWDKQIQSLCFQVNNLLEKISQTAPEWT'
        'AQAMEAQMAQ'},
    'Q92905': {'db_code':'CSN5_HUMAN', 'sequence':
        'MAASGSGMAQKTWELANNMQEAQSIDEIYKYDKKQQQEILAAKPWTKDHHYFKYCKISALALLKMV'
        'MHARSGGNLEVMGLMLGKVDGETMIIMDSFALPVEGTETRVNAQAAAYEYMAAYIENAKQVGRLEN'
        'AIGWYHSHPGYGCWLSGIDVSTQMLNQQFQEPFVAVVIDPTRTISAGKVNLGAFRTYPKGYKPPDE'
        'GPSEYQTIPLNKIEDFGVHCKQYYALEVSYFKSSLDRKLLELLWNKYWVNTLSSSSLLTNADYTTG'
        'QVFDLSEKLEQSEAQLGRGSFMLGLETHDRKSEDKLAKATRDSCKTTIEAIHGLMSQVIKDKLFNQ'
        'INIS'},
    'Q7L5N1': {'db_code':'CSN6_HUMAN', 'sequence':
        'MAAAAAAAAATNGTGGSSGMEVDAAVVPSVMACGVTGSVSVALHPLVILNISDHWIRMRSQEGRPV'
        'QVIGALIGKQEGRNIEVMNSFELLSHTVEEKIIIDKEYYYTKEEQFKQVFKELEFLGWYTTGGPPD'
        'PSDIHVHKQVCEIIESPLFLKLNPMTKHTDLPVSVFESVIDIINGEATMLFAELTYTLATEEAERI'
        'GVDHVARMTATGSGENSTVAEHLIAQHSAIKMLHSRVKLILEYVKASEAGEVPFNHEILREAYALC'
        'HCLPVLSTDKFKTDFYDQCNDVGLMAYLGTITKTCNTMNQFVNKFNVLYDRQGIGRRMRGLFF'},
    'Q9H9Q2': {'db_code':'CSN7B_HUMAN', 'sequence':
        'MAGEQKPSSNLLEQFILLAKGTSGSALTALISQVLEAPGVYVFGELLELANVQELAEGANAAYLQL'
        'LNLFAYGTYPDYIANKESLPELSTAQQNKLKHLTIVSLASRMKCIPYSVLLKDLEMRNLRELEDLI'
        'IEAVYTDIIQGKLDQRNQLLEVDFCIGRDIRKKDINNIVKTLHEWCDGCEAVLLGIEQQVLRANQY'
        'KENHNRTQQQVEAEVTNIKKTLKATASSSAQEMEQQLAERECPPHAEQRQPTKKMSKVKGLVSSRH'},
    'Q99627': {'db_code':'CSN8_HUMAN', 'sequence':
        'MPVAVMAESAFSFKKLLDQCENQELEAPGGIATPPVYGQLLALYLLHNDMNNARYLWKRIPPAIKS'
        'ANSELGGIWSVGQRIWQRDFPGIYTTINAHQWSETVQPIMEALRDATRRRAFALVSQAYTSIIADD'
        'FAAFVGLPVEEAVKGILEQGWQADSTTRMVLPRKPVAGALDVSFNKFIPLSEPAPVPPIPNEQQLA'
        'RLTDYVAFLEN'},
    'Q8WXC6': {'db_code':'CSN9_HUMAN', 'sequence':
        'MKPAVDEMFPEGAGPYVDLDEAGGSTGLLMDLAANEKAVHADFFNDFEDLFDDDDIQ'},
}

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
    if hasattr(ihm, 'reference'):
        def from_acc(cls, accession):
            return cls(accession=accession, **accession_map[accession])
        ihm.reference.UniProtSequence.from_accession = classmethod(from_acc)
