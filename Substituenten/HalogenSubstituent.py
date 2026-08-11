from InputParser import dic_converter
from HelferGeometrie import _substituent_verbindung

def halogen_substituent(stamm_kette_punkte, halogen_input, plotter, bindung_verschiebung, besetzt_liste):
    try:
        # Name, Farbe, Grösse
        halogen_zeichnung = {
            "fluor": ('F', "#FFD1DC", 40),
            "chlor": ('Cl', "#228B22", 50),
            "brom": ('Br', "#CC5500", 60),
            "iod": ('I', "#9D00FF", 70)
        }

        alle_halogene_alle_pos = dic_converter(halogen_input)

        if alle_halogene_alle_pos is None:
            alle_halogene_alle_pos = {}

        endpunkt_liste = _substituent_verbindung(stamm_kette_punkte, alle_halogene_alle_pos, plotter,
                                                bindung_verschiebung, besetzt_liste)

        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            halogen = endpunkt[2]
            halogen_lower = halogen.lower()

            text, color, point_size = halogen_zeichnung[halogen_lower]

            plotter.add_point_labels(
                points=[endpunkt[0]],
                labels=[text],
                font_size=40,
                point_color=color,
                point_size=point_size,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )
    except Exception as e:
        print(f"Fehler in der Darstellung der Elemente:{e}")
