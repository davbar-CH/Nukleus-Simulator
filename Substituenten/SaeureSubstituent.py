from HelferGeometrie import _saeure_endpunkt_zeichnen
from InputParser import dic_converter

def saeure_substituent(stamm_kette_punkte, saeure_input, plotter, bindung_verschiebung, besetzt_liste, verschiebung_h=0.2):
    try:
        alle_saeure_alle_pos = dic_converter(saeure_input)
        zaehl_wort = saeure_input[0][1]
        letzte_position = len(stamm_kette_punkte)
        verschiebung_sauerstoff = bindung_verschiebung
        if zaehl_wort == "di":
            alle_saeure_alle_pos["säure"] = [letzte_position]
            verschiebung_alkohol_ende = (bindung_verschiebung + 90 if letzte_position % 2 == 0
                                  else bindung_verschiebung - 180)
            _saeure_endpunkt_zeichnen(stamm_kette_punkte, letzte_position, alle_saeure_alle_pos,
                                       plotter, verschiebung_sauerstoff, verschiebung_alkohol_ende, besetzt_liste)

        alle_saeure_alle_pos["säure"] = [1]
        verschiebung_alkohol_ende = bindung_verschiebung - 90
        _saeure_endpunkt_zeichnen(stamm_kette_punkte, 1, alle_saeure_alle_pos,
                                   plotter, verschiebung_sauerstoff, verschiebung_alkohol_ende, besetzt_liste)

    except Exception as e:
        print(f"Fehler in der Darstellung der Säure / des Aldehyds: {e}")
