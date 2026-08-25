from HelferGeometrie import *


def keton_substituent(stamm_kette_punkte, alle_saeure_alle_pos,
                      plotter, verschiebung_sauerstoff,
                      besetzt_liste):

    endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_saeure_alle_pos, plotter,
                                          verschiebung_sauerstoff, besetzt_liste)
    for endpunkt in endpunkt_liste:
        sauerstoff_koordinaten = endpunkt[0]
        sub_pos = endpunkt[1]

        verbindung_koordinaten = [sauerstoff_koordinaten, stamm_kette_punkte[sub_pos]]

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
