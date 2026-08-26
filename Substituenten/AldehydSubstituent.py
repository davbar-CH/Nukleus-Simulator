from InputParser import dic_converter
from Substituenten import keton_substituent
from HelferGeometrie import substituent_verbindung

def aldehyd_substituent(stamm_kette_punkte, aldehyd_input, plotter, bindung_verschiebung, besetzt_liste, verschiebung_h=0.2):
    try:
        alle_aldehyd_alle_pos = dic_converter(aldehyd_input)
        zaehl_wort = aldehyd_input[0][1]
        letzte_position = len(stamm_kette_punkte)
        verschiebung_sauerstoff = bindung_verschiebung

        if zaehl_wort == "di":
            alle_aldehyd_alle_pos["al"] = [letzte_position]
            verschiebung_wasserstoff_ende = (bindung_verschiebung + 90 if letzte_position % 2 == 0
                                  else bindung_verschiebung - 180)

            endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_aldehyd_alle_pos, plotter,
                                                    verschiebung_wasserstoff_ende, besetzt_liste)
            if endpunkt_liste is None:
                endpunkt_liste = []

            for endpunkt in endpunkt_liste:
                wasserstoff_anfangspunkt = endpunkt[0]

                plotter.add_point_labels(
                    points=wasserstoff_anfangspunkt,
                    labels=["H"],
                    font_size=40,
                    point_color="#d9e4ea",
                    point_size=20,
                    render_points_as_spheres=True,
                    always_visible=True,
                    shape=None
                )

            besetzt_liste.remove(letzte_position)
            keton_substituent(stamm_kette_punkte, alle_aldehyd_alle_pos,
                                       plotter, verschiebung_sauerstoff, besetzt_liste)

        position = 1
        if "formyl" in alle_aldehyd_alle_pos:
            alle_aldehyd_alle_pos["formyl"] = [letzte_position]
            position = letzte_position
            verschiebung_wasserstoff_anfang = (bindung_verschiebung + 90 if letzte_position % 2 == 0
                                  else bindung_verschiebung - 180)
        else:
            alle_aldehyd_alle_pos["al"] = [position]
            verschiebung_wasserstoff_anfang = bindung_verschiebung + 180

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_aldehyd_alle_pos, plotter,
                                                verschiebung_wasserstoff_anfang, besetzt_liste)
        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            wasserstoff_anfangspunkt = endpunkt[0]

            plotter.add_point_labels(
                points=wasserstoff_anfangspunkt,
                labels=["H"],
                font_size=40,
                point_color="#d9e4ea",
                point_size=20,
                render_points_as_spheres=True,
                always_visible=True,
                shape=None
            )
        besetzt_liste.remove(position)
        keton_substituent(stamm_kette_punkte, alle_aldehyd_alle_pos,
                                   plotter, verschiebung_sauerstoff, besetzt_liste)

    except Exception as e:
        print(f"Fehler in der Darstellung der Aldehyds: {e}")