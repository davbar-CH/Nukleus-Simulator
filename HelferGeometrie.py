from pyvista import lines_from_points
from numpy import array, sin, cos, linalg


def substituent_verbindung(stamm_kette_punkte, substituent_dic, plotter, bindung_verschiebung, besetzt_liste):
    try:
        endpunkt_liste = []

        if substituent_dic is not None:
            for substituent in substituent_dic:
                for sub_pos in substituent_dic.get(substituent, {}):
                    anfangspunkt = stamm_kette_punkte[sub_pos - 1]
                    vorzeichen = -1 if sub_pos in besetzt_liste else 1
                    y_formel = -1 if sub_pos % 2 == 0 else 1

                    substituent_verbindung_punkte = array([
                        anfangspunkt,
                        array([anfangspunkt[0] + vorzeichen * (0.5 * cos(bindung_verschiebung)),
                               anfangspunkt[1] + vorzeichen * y_formel * (0.5 * sin(bindung_verschiebung)),
                               0])
                    ])

                    endpunkt = [substituent_verbindung_punkte[1], sub_pos - 1, substituent]
                    endpunkt_liste.append(endpunkt)

                    besetzt_liste.append(sub_pos)
                    verbindung_substituent = lines_from_points(substituent_verbindung_punkte)
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

