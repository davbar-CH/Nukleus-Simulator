import re

def dic_converter(input):
    """
    :param input: Der jeweilige Text in der Form ('2,3', 'Di', 'ethyl'), ('4,5', 'di', 'methyl')
    :return: {'ethyl': [2, 3], 'methyl': [5, 5]}, Zahlen sind danach integer
    """

    try:
        alle_pos = [x[0] for x in input]
        alle_gruppen = [x[2].lower() for x in input]

        alle_gruppen_alle_pos = {}
        for i, position in enumerate(alle_pos):
            bindung_pos = {alle_gruppen[i]: [int(x) for x in re.findall(r"\d", position)]}
            alle_gruppen_alle_pos.update(bindung_pos)

        return alle_gruppen_alle_pos
    except Exception as e:
        print(f"Fehler in der dictionary Konvertierung:{e}")

def text_auslesen(input_text):
    """
    Die Funktion liest mit Regex den Text aus und gibt alle Substituenten, Bindungen und Stereoisomerie zurück

    :param input_text: Text aus der Textbox, in der Form "(Stereo)-Position-Substituent(en)-Stamm-Position-Bindung".
    :return: gibt alle Substituenten, Bindungen und Stereoisomerie zurück.
    """
    zaehl_woerter_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?"

    stereo_pattern = r"\(([EZRSezrs,\d]+)\)-"

    pattern_halogen = rf"{zaehl_woerter_pattern}(fluor|chlor|brom|iod)"
    pattern_alkan = rf"{zaehl_woerter_pattern}(methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|octyl|nonyl|decyl)"
    pattern_phenyl = rf"{zaehl_woerter_pattern}(phenyl)"
    amin_pattern = rf"{zaehl_woerter_pattern}(amin|amino)"
    alkohol_pattern = rf"{zaehl_woerter_pattern}(ol|hydroxy)"

    alle_sub_pattern = rf"{zaehl_woerter_pattern}(fluor|chlor|brom|iod|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl\
    |phenyl|hydroxy|amino)"

    stereo = re.findall(stereo_pattern, input_text, flags=re.IGNORECASE)

    halogen = re.findall(pattern_halogen, input_text, flags=re.IGNORECASE)
    alkan = re.findall(pattern_alkan, input_text, flags=re.IGNORECASE)
    phenyl = re.findall(pattern_phenyl, input_text, flags=re.IGNORECASE)
    amin = re.findall(amin_pattern, input_text, flags=re.IGNORECASE)
    alkohol = re.findall(alkohol_pattern, input_text, flags=re.IGNORECASE)

    input_ohne_stereo = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)
    input_ohne_stereo_sub = re.sub(alle_sub_pattern, '', input_ohne_stereo, flags=re.IGNORECASE)

    is_cyclo = False
    if re.search(r"cyclo", input_ohne_stereo_sub, re.IGNORECASE):
        is_cyclo = True

    stamm_pattern = r"(?:cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
    stamm = re.findall(stamm_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    bindung_pattern = rf"{zaehl_woerter_pattern}(en|in)(?!\w)"
    bindung_typ = re.findall(bindung_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    saeure_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di))?(säure)"
    saeure = re.findall(saeure_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    aldehyd_pattern = rf"{zaehl_woerter_pattern}(formyl|al)"
    aldehyd = re.findall(aldehyd_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    keton_pattern = rf"{zaehl_woerter_pattern}(oxo|on)"
    keton = re.findall(keton_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    return stereo, alkan, halogen, phenyl, alkohol, amin, is_cyclo, stamm, bindung_typ, saeure, aldehyd, keton
