from pyvista import lines_from_points
from numpy import array
from InputParser import dic_converter
from HelferGeometrie import substituent_verbindung

def amino_substituent(stamm_kette_punkte, amin_input, plotter,
                      bindung_verschiebung, besetzt_liste, verschiebung_h=0.2):
    try:
        alle_amin_alle_pos = dic_converter(amin_input)

        if alle_amin_alle_pos is None:
            alle_amin_alle_pos = {}

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_amin_alle_pos, plotter, bindung_verschiebung, besetzt_liste)

        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            wasserstoff_anfangspunkt = endpunkt[0]
            position_kette = endpunkt[1]

            plotter.add_point_labels(
                points=wasserstoff_anfangspunkt,
                labels=["N"],
                font_size=30,
                point_color="#1e7fcb",
                point_size=40,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )

            wasserstoff_endpunkt = 2 * wasserstoff_anfangspunkt - stamm_kette_punkte[position_kette]

            wasserstoff_endpunkt_links = [wasserstoff_endpunkt[0] + verschiebung_h, wasserstoff_endpunkt[1], wasserstoff_endpunkt[2]]
            wasserstoff_endpunkt_rechts = [wasserstoff_endpunkt[0] - verschiebung_h, wasserstoff_endpunkt[1], wasserstoff_endpunkt[2]]


            wasserstoff_verbindung_linie_links = lines_from_points([wasserstoff_anfangspunkt, wasserstoff_endpunkt_links])
            wasserstoff_verbindung_linie_rechts = lines_from_points([wasserstoff_anfangspunkt, wasserstoff_endpunkt_rechts])

            plotter.add_mesh(wasserstoff_verbindung_linie_links, line_width=3)
            plotter.add_mesh(wasserstoff_verbindung_linie_rechts, line_width=3)

            plotter.add_point_labels(
                points=[wasserstoff_endpunkt_links, wasserstoff_endpunkt_rechts],
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
