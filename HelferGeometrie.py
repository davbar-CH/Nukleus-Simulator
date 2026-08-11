from pyvista import lines_from_points
from numpy import array, sin, cos, linalg
from Substituenten.AlkoholSubstituent import alkohol_substituent

def _substituent_verbindung(stamm_kette_punkte, substituent_dic, plotter, bindung_verschiebung, besetzt_liste):
    try:
        endpunkt_liste = []

        if substituent_dic is not None:
            for substituent in substituent_dic:
                for sub_pos in substituent_dic.get(substituent, {}):
                    anfangspunkt = stamm_kette_punkte[sub_pos - 1]
                    vorzeichen = -1 if sub_pos in besetzt_liste else 1
                    y_formel = -1 if sub_pos % 2 == 0 else 1

                    _substituent_verbindung_punkte = array([
                        anfangspunkt,
                        array([anfangspunkt[0] + vorzeichen * (0.5 * cos(bindung_verschiebung)),
                                  anfangspunkt[1] + vorzeichen * y_formel * (0.5 * sin(bindung_verschiebung)),
                                  0])
                    ])

                    endpunkt = [_substituent_verbindung_punkte[1], sub_pos - 1, substituent]
                    endpunkt_liste.append(endpunkt)

                    besetzt_liste.append(sub_pos)
                    verbindung_substituent = lines_from_points(_substituent_verbindung_punkte)
                    plotter.add_mesh(verbindung_substituent, line_width=2, color=(0, 0, 0))
            return endpunkt_liste

        else:
            return []

    except Exception as e:
        print(f"Fehler in der Darstellung der Substituent-Verbindung: {e}")

def bindung_zeichnen(stamm_kette_punkte, plotter, alle_bindungen_alle_pos, verschiebung_bindung=0.1,
                     laenge_bindung=0.2):
    try:
        for bindung in alle_bindungen_alle_pos:
            for bindung_pos in alle_bindungen_alle_pos.get(bindung):
                p1 = array(stamm_kette_punkte[bindung_pos - 1][:2])
                p2 = array(stamm_kette_punkte[bindung_pos][:2])

                richtung = p2 - p1

                normale = array([-richtung[1], richtung[0]])
                normale = normale / linalg.norm(normale)

                p1_verschoben_oben = (p1 + verschiebung_bindung * normale) + laenge_bindung * richtung
                p2_verschoben_oben = (p2 + verschiebung_bindung * normale) - laenge_bindung * richtung

                alken_punkte = array([
                    array([p1_verschoben_oben[0], p1_verschoben_oben[1], 0]),
                    array([p2_verschoben_oben[0], p2_verschoben_oben[1], 0])
                ])
                if bindung == "in":
                    p1_verschoben_unten = (p1 - verschiebung_bindung * normale) + laenge_bindung * richtung
                    p2_verschoben_unten = (p2 - verschiebung_bindung * normale) - laenge_bindung * richtung

                    alkin_punkte = array([
                        array([p1_verschoben_unten[0], p1_verschoben_unten[1], 0]),
                        array([p2_verschoben_unten[0], p2_verschoben_unten[1], 0])
                    ])

                    alkin_kette = lines_from_points(alkin_punkte)
                    plotter.add_mesh(alkin_kette, line_width=2, color=(255, 0, 0))

                alken_kette = lines_from_points(alken_punkte)
                plotter.add_mesh(alken_kette, line_width=2, color=(255, 0, 0))

    except Exception as e:
        print(f"Fehler in der Darstellung der Bindung:{e}")

def _saeure_endpunkt_zeichnen(stamm_kette_punkte, position, alle_saeure_alle_pos,
                              plotter, verschiebung_sauerstoff, verschiebung_alkohol_ende,
                              besetzt_liste):
    """Setzt den OH-Substituenten an `position`, verbindet ihn mit dem Stamm
    und zeichnet die C=O-Doppelbindung samt O-Label."""
    alkohol_substituent(stamm_kette_punkte, [(str(position), "", "ol")], plotter,
                         verschiebung_alkohol_ende, besetzt_liste)
    besetzt_liste.remove(position)

    endpunkt_liste = _substituent_verbindung(stamm_kette_punkte, alle_saeure_alle_pos, plotter,
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
