from InputParser import dic_converter
from HelferGeometrie import substituent_verbindung, bindung_zeichnen
from pyvista import Polygon

def phenyl_substituent(stamm_kette_punkte, phenyl_input, plotter, bindung_verschiebung, besetzt_liste, phenyl_groesse=0.5):
    try:
        alle_phenyl_sub_alle_pos = dic_converter(phenyl_input)

        if alle_phenyl_sub_alle_pos is None:
            alle_phenyl_sub_alle_pos = {}

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_phenyl_sub_alle_pos, plotter,
                                                bindung_verschiebung, besetzt_liste)
        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            koordinaten = endpunkt[0]
            position_kette = endpunkt[1]

            phenyl_mesh = Polygon(
                center=[koordinaten[0],
                        koordinaten[1] - phenyl_groesse if position_kette % 2 == 0 else koordinaten[1] + phenyl_groesse,
                        0],
                radius=phenyl_groesse,
                fill=False)
            plotter.add_mesh(phenyl_mesh, line_width=2, color=(0, 0, 0))
            punkte = phenyl_mesh.points

            bindung_zeichnen(punkte, plotter, alle_bindungen_alle_pos={"en": (1, 3, 5)}, verschiebung_bindung=-0.05,
                             laenge_bindung=0.008)

    except Exception as e:
        print(f"Fehler in der Darstellung von Phenyl: {e}")
