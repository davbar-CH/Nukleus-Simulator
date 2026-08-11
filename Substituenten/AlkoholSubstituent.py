from pyvista import lines_from_points
from InputParser import dic_converter
from HelferGeometrie import _substituent_verbindung

def alkohol_substituent(stamm_kette_punkte, alkohol_input, plotter,
                        bindung_verschiebung, besetzt_liste, verschiebung_h=0.2):
    try:
        alle_alkohol_alle_pos = dic_converter(alkohol_input)

        if alle_alkohol_alle_pos is None:
            alle_alkohol_alle_pos = {}

        endpunkt_liste = _substituent_verbindung(stamm_kette_punkte, alle_alkohol_alle_pos, plotter,
                                                bindung_verschiebung, besetzt_liste)
        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            wasserstoff_anfangspunkt = endpunkt[0]
            position_kette = endpunkt[1]

            plotter.add_point_labels(
                points=wasserstoff_anfangspunkt,
                labels=["O"],
                font_size=40,
                point_color="#ec0c0d",
                point_size=40,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )

            wasserstoff_endpunkt = 2 * wasserstoff_anfangspunkt - stamm_kette_punkte[position_kette]

            wasserstoff_verbindung_linie = lines_from_points([wasserstoff_anfangspunkt, wasserstoff_endpunkt])
            plotter.add_mesh(wasserstoff_verbindung_linie, line_width=3)

            plotter.add_point_labels(
                points=[wasserstoff_endpunkt],
                labels=["H"],
                font_size=40,
                point_color="#d9e4ea",
                point_size=20,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )

    except Exception as e:
        print(f"Fehler in der Darstellung vom Alkohol: {e}")
