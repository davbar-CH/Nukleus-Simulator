from pyvista import lines_from_points
from numpy import array
from InputParser import dic_converter
from HelferGeometrie import _substituent_verbindung

def amino_substituent(stamm_kette_punkte, amin_input, plotter,
                      bindung_verschiebung, besetzt_liste, verschiebung_h=0.2):
    try:
        alle_amin_alle_pos = dic_converter(amin_input)

        if alle_amin_alle_pos is None:
            alle_amin_alle_pos = {}

        endpunkt_liste = _substituent_verbindung(stamm_kette_punkte, alle_amin_alle_pos, plotter, bindung_verschiebung, besetzt_liste)

        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            koordinaten = endpunkt[0]
            position_kette = endpunkt[1]

            plotter.add_point_labels(
                points=koordinaten,
                labels=["N"],
                font_size=30,
                point_color="#1e7fcb",
                point_size=40,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )

            vorzeichen = -1 if position_kette in besetzt_liste else 1
            y_formel = bindung_verschiebung - 2 if position_kette % 2 == 0 else bindung_verschiebung

            wasserstoff_verbindung_punkte_links = array([
                koordinaten,
                array([stamm_kette_punkte[position_kette][0] + vorzeichen * -bindung_verschiebung,
                          y_formel - verschiebung_h if position_kette % 2 == 0 else y_formel + verschiebung_h, 0])
            ])

            wasserstoff_verbindung_punkte_rechts = array([
                koordinaten,
                array([2 * koordinaten[0] - wasserstoff_verbindung_punkte_links[1][0],
                          y_formel - verschiebung_h if position_kette % 2 == 0 else y_formel + verschiebung_h, 0])
            ])

            wasserstoff_verbindung_linie_links = lines_from_points(wasserstoff_verbindung_punkte_links)
            wasserstoff_verbindung_linie_rechts = lines_from_points(wasserstoff_verbindung_punkte_rechts)

            plotter.add_mesh(wasserstoff_verbindung_linie_links, line_width=3)
            plotter.add_mesh(wasserstoff_verbindung_linie_rechts, line_width=3)

            plotter.add_point_labels(
                points=[wasserstoff_verbindung_punkte_links[1], wasserstoff_verbindung_punkte_rechts[1]],
                labels=["H", "H"],
                font_size=30,
                point_color="#d9e4ea",
                point_size=20,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )

    except Exception as e:
        print(f"Fehler in der Darstellung vom Amin: {e}")
