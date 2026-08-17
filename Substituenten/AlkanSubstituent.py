from pyvista import lines_from_points
from numpy import array, insert
from InputParser import dic_converter


def alkan_substituent(stamm_kette_punkte, alkan_input, plotter, bindung_verschiebung, besetzt_liste, dehnung_x=8, dehnung_y=0.5):
    try:
        """
        pos = Position des Substituents, in der Form 1 oder 3,8
        sub = Substituent, in der Form Methyl oder methyl

        Ist eine Position besetzt, dann wird das Vorzeichen gekehrt, 
        damit die Substituenten trotzdem angezeigt werden können

        Bei ungeraden Positionen wird bei y=0 gestartet, bei geraden Positionen bei y=1

        Der Substituent startet bei der jeweiligen x-Koordinate (pos-1). Die x-Koordinate alterniert jeweils zwischen 
        der x-Koordinate und x-Koordinate - 0.25. 
        Die y-Koordinate wird stets um 0.5 grösser (bei geraden Positionen) bzw.
        um 0.5 kleiner bei ungeraden Positionen.
        """
        substituent_laenge = {
            "methyl": 1,
            "ethyl": 2,
            "propyl": 3,
            "butyl": 4,
            "pentyl": 5,
            "hexyl": 6,
            "heptyl": 7,
            "octyl": 8,
            "nonyl": 9,
            "decyl": 10,
        }

        alle_alkane_alle_pos = dic_converter(alkan_input)

        if alle_alkane_alle_pos is None:
            alle_alkane_alle_pos = {}

        for alkan in alle_alkane_alle_pos:
            for sub_pos in alle_alkane_alle_pos.get(alkan, {}):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = (lambda x: (x + 2) * dehnung_y) if sub_pos % 2 == 0 else (lambda x: x * -dehnung_y)
                laenge = substituent_laenge.get(alkan) + 1

                anfangspunkt = stamm_kette_punkte[sub_pos - 1]

                sub_alkan_punkte = array([
                    [stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * ((-1 + (-1) ** x) / dehnung_x), y_formel(x), 0]
                    for x in range(1, laenge)
                ])

                sub_alkan_punkte = insert(sub_alkan_punkte, 0, anfangspunkt, axis=0)
                besetzt_liste.append(sub_pos)
                alkan_kette = lines_from_points(sub_alkan_punkte)
                plotter.add_mesh(alkan_kette, line_width=2, color=(0, 0, 0))

    except Exception as e:
        print(f"Kein Alkan-Substituent: {e}")
