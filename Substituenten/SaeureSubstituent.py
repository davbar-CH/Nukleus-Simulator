from HelferGeometrie import *
from InputParser import dic_converter
from Substituenten import alkohol_substituent

def sauerstoff_doppelbindung(stamm_kette_punkte, position, alle_saeure_alle_pos,
                              plotter, verschiebung_sauerstoff,
                              besetzt_liste):
    besetzt_liste.remove(position)

    endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_saeure_alle_pos, plotter,
                                             verschiebung_sauerstoff, besetzt_liste)
    sauerstoff_koordinaten = endpunkt_liste[0][0]
    verbindung_koordinaten = [sauerstoff_koordinaten, stamm_kette_punkte[position - 1]]

    bindung_zeichnen(verbindung_koordinaten, plotter, alle_bindungen_alle_pos={"en": [1]},
                     verschiebung_bindung=-0.05, laenge_bindung=0.008)

    plotter.add_point_labels(
        points=sauerstoff_koordinaten,
        labels=["O"],
        font_size=40,
        point_color="#ec0c0d",
        point_size=40,
        render_points_as_spheres=True,
        always_visible=True,
        shape=None
    )

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
            alkohol_substituent(stamm_kette_punkte, [(str(letzte_position), "", "ol")], plotter,
                                verschiebung_alkohol_ende, besetzt_liste)
            sauerstoff_doppelbindung(stamm_kette_punkte, letzte_position, alle_saeure_alle_pos,
                                       plotter, verschiebung_sauerstoff, besetzt_liste)

        alle_saeure_alle_pos["säure"] = [1]
        verschiebung_alkohol_ende = bindung_verschiebung - 90
        alkohol_substituent(stamm_kette_punkte, [(str(1), "", "ol")], plotter,
                            verschiebung_alkohol_ende, besetzt_liste)
        sauerstoff_doppelbindung(stamm_kette_punkte, 1, alle_saeure_alle_pos,
                                   plotter, verschiebung_sauerstoff, besetzt_liste)

    except Exception as e:
        print(f"Fehler in der Darstellung der Säure: {e}")
