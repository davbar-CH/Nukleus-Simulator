import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import *
import numpy as np
import re


def text_auslesen(input_text):
    """
    Die Funktion liest mit Regex den Text aus und gibt alle Substituenten, Bindungen und Stereoisomerie zurück

    :param input_text: Text aus der Textbox, in der Form "(Stereo)-Position-Substituent(en)-Stamm-Position-Bindung".
    :return: gibt alle Substituenten, Bindungen und Stereoisomerie zurück.
    """
    zaehl_woerter_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?"

    stereo_pattern = r"\(([EZRSezrs,\d]+)\)-"

    pattern_halogen = rf"{zaehl_woerter_pattern}(fluor|chlor|brom|iod)"
    pattern_alkan = rf"{zaehl_woerter_pattern}(methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|octyl|nonyl|decyl)"
    pattern_phenyl = rf"{zaehl_woerter_pattern}(phenyl)"
    amin_pattern = rf"{zaehl_woerter_pattern}(amin|amino)"
    alkohol_pattern = rf"{zaehl_woerter_pattern}(ol|hydroxy)"

    alle_sub_pattern = rf"{zaehl_woerter_pattern}(fluor|chlor|brom|iod|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl\
    |phenyl|hydroxy|amino)"

    stereo = re.findall(stereo_pattern, input_text, flags=re.IGNORECASE)

    halogen = re.findall(pattern_halogen, input_text, flags=re.IGNORECASE)
    alkan = re.findall(pattern_alkan, input_text, flags=re.IGNORECASE)
    phenyl = re.findall(pattern_phenyl, input_text, flags=re.IGNORECASE)
    amin = re.findall(amin_pattern, input_text, flags=re.IGNORECASE)
    alkohol = re.findall(alkohol_pattern, input_text, flags=re.IGNORECASE)

    input_ohne_stereo = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)
    input_ohne_stereo_sub = re.sub(alle_sub_pattern, '', input_ohne_stereo, flags=re.IGNORECASE)

    is_cyclo = False
    if re.search(r"cyclo", input_ohne_stereo_sub, re.IGNORECASE):
        is_cyclo = True

    stamm_pattern = r"(?:cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
    stamm = re.findall(stamm_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    bindung_pattern = rf"{zaehl_woerter_pattern}(en|in)(?!\w)"
    bindung_typ = re.findall(bindung_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    saeure_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di))?(säure)"
    saeure = re.findall(saeure_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    aldehyd_pattern = rf"{zaehl_woerter_pattern}(formyl|al)"
    aldehyd = re.findall(aldehyd_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    return stereo, alkan, halogen, phenyl, alkohol, amin, is_cyclo, stamm, bindung_typ, saeure, aldehyd


def dic_converter(input):
    """
    :param input: Der jeweilige Text in der Form ('2,3', 'Di', 'ethyl'), ('4,5', 'di', 'methyl')
    :return: {'ethyl': [2, 3], 'methyl': [5, 5]}, Zahlen sind danach integer
    """

    try:
        alle_pos = [x[0] for x in input]
        alle_gruppen = [x[2].lower() for x in input]

        alle_gruppen_alle_pos = {}
        for i, position in enumerate(alle_pos):
            bindung_pos = {alle_gruppen[i]: [int(x) for x in re.findall(r"\d", position)]}
            alle_gruppen_alle_pos.update(bindung_pos)

        return alle_gruppen_alle_pos
    except Exception as e:
        print(f"Fehler in der dictionary Konvertierung:{e}")


def stamm_kette(stamm_input, dehnung_x=0.5, dehnung_y=2):
    try:
        stamm_laenge = {
            "eth": 2,
            "prop": 3,
            "but": 4,
            "pent": 5,
            "hex": 6,
            "hept": 7,
            "oct": 8,
            "non": 9,
            "dec": 10,
        }

        global besetzt_liste
        besetzt_liste = []
        stamm_kette_punkte = []

        stamm = stamm_input[0].lower()
        laenge = stamm_laenge.get(stamm)
        if laenge is not None:
            stamm_kette_punkte = np.array(
                [[x * dehnung_x, (1 - (-1) ** x) / dehnung_y, 0] for x in range(0, laenge)])
        else:
            print(f"Keine solche Stammkette{stamm}")

        return stamm_kette_punkte

    except Exception as e:
        print(f"Fehler in der stamm_kette: {e}")


def substituent_verbindung(stamm_kette_punkte, substituent_dic, plotter, bindung_verschiebung):
    try:
        endpunkt_liste = []

        if substituent_dic is not None:
            for substituent in substituent_dic:
                for sub_pos in substituent_dic.get(substituent, {}):
                    anfangspunkt = stamm_kette_punkte[sub_pos - 1]
                    vorzeichen = -1 if sub_pos in besetzt_liste else 1
                    y_formel = -1 if sub_pos % 2 == 0 else 1

                    substituent_verbindung_punkte = np.array([
                        anfangspunkt,
                        np.array([anfangspunkt[0] + vorzeichen * (0.5 * np.cos(bindung_verschiebung)),
                                  anfangspunkt[1] + vorzeichen * y_formel * (0.5 * np.sin(bindung_verschiebung)),
                                  0])
                    ])

                    endpunkt = [substituent_verbindung_punkte[1], sub_pos - 1, substituent]
                    endpunkt_liste.append(endpunkt)

                    besetzt_liste.append(sub_pos)
                    verbindung_substituent = pv.lines_from_points(substituent_verbindung_punkte)
                    plotter.add_mesh(verbindung_substituent, line_width=2, color=(0, 0, 0))
            return endpunkt_liste

        else:
            return []


    except Exception as e:
        print(f"Fehler in der Darstellung der Substituent-Verbindung: {e}")


def stamm_anpassung(stamm_kette_punkte, plotter, alle_bindungen_alle_pos):
    try:
        alle_positionen_en = alle_bindungen_alle_pos.get("en", {})
        alle_positionen_in = alle_bindungen_alle_pos.get("in", {})

        if alle_positionen_en:
            for i, pos in enumerate(alle_positionen_en):

                if i + 1 == len(alle_positionen_en):
                    break

                if alle_positionen_en[i + 1] - pos == 1:
                    if stamm_kette_punkte[pos][1] == 0:
                        stamm_kette_punkte[pos][1] = 1
                    elif stamm_kette_punkte[pos][1] == 1:
                        stamm_kette_punkte[pos][1] = 0
                else:
                    pass

        elif alle_positionen_in:
            for pos in alle_positionen_in:

                if pos - 1 == 0:
                    new_val = 0

                else:
                    new_val = 1 - stamm_kette_punkte[pos - 1][1]

                    for j in range(pos + 1, len(stamm_kette_punkte)):
                        stamm_kette_punkte[j][1] = new_val if j == pos + 1 else 1 - stamm_kette_punkte[j][1]

                stamm_kette_punkte[pos - 1][1] = new_val
                stamm_kette_punkte[pos][1] = new_val

        stamm_kette = pv.lines_from_points(stamm_kette_punkte)
        plotter.add_mesh(stamm_kette, line_width=4, color=(0, 0, 0))

        return alle_bindungen_alle_pos

    except Exception as e:
        print(f"Fehler in der stamm_anpassung: {e}")


def bindung_zeichnen(stamm_kette_punkte, plotter, alle_bindungen_alle_pos, verschiebung_bindung=0.1,
                     laenge_bindung=0.2):
    try:
        for bindung in alle_bindungen_alle_pos:
            for bindung_pos in alle_bindungen_alle_pos.get(bindung):
                p1 = np.array(stamm_kette_punkte[bindung_pos - 1][:2])
                p2 = np.array(stamm_kette_punkte[bindung_pos][:2])

                richtung = p2 - p1

                normale = np.array([-richtung[1], richtung[0]])
                normale = normale / np.linalg.norm(normale)

                p1_verschoben_oben = (p1 + verschiebung_bindung * normale) + laenge_bindung * richtung
                p2_verschoben_oben = (p2 + verschiebung_bindung * normale) - laenge_bindung * richtung

                alken_punkte = np.array([
                    np.array([p1_verschoben_oben[0], p1_verschoben_oben[1], 0]),
                    np.array([p2_verschoben_oben[0], p2_verschoben_oben[1], 0])
                ])
                if bindung == "in":
                    p1_verschoben_unten = (p1 - verschiebung_bindung * normale) + laenge_bindung * richtung
                    p2_verschoben_unten = (p2 - verschiebung_bindung * normale) - laenge_bindung * richtung

                    alkin_punkte = np.array([
                        np.array([p1_verschoben_unten[0], p1_verschoben_unten[1], 0]),
                        np.array([p2_verschoben_unten[0], p2_verschoben_unten[1], 0])
                    ])

                    alkin_kette = pv.lines_from_points(alkin_punkte)
                    plotter.add_mesh(alkin_kette, line_width=2, color=(255, 0, 0))

                alken_kette = pv.lines_from_points(alken_punkte)
                plotter.add_mesh(alken_kette, line_width=2, color=(255, 0, 0))

    except Exception as e:
        print(f"Fehler in der Darstellung der Bindung:{e}")


def alkan_substituent(stamm_kette_punkte, alkan_input, plotter, dehnung_x=8, dehnung_y=0.5):
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

                sub_alkan_punkte = np.array([
                    [stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * ((-1 + (-1) ** x) / dehnung_x), y_formel(x), 0]
                    for x in range(1, laenge)
                ])

                sub_alkan_punkte = np.insert(sub_alkan_punkte, 0, anfangspunkt, axis=0)
                besetzt_liste.append(sub_pos)
                alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                plotter.add_mesh(alkan_kette, line_width=2, color=(0, 0, 0))

    except Exception as e:
        print(f"Kein Alkan-Substituent: {e}")


def halogen_substituent(stamm_kette_punkte, halogen_input, plotter, bindung_verschiebung):
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

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_halogene_alle_pos, plotter,
                                                bindung_verschiebung)

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


def phenyl_substituent(stamm_kette_punkte, phenyl_input, plotter, bindung_verschiebung, phenyl_groesse=0.5):
    try:
        alle_phenyl_sub_alle_pos = dic_converter(phenyl_input)

        if alle_phenyl_sub_alle_pos is None:
            alle_phenyl_sub_alle_pos = {}

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_phenyl_sub_alle_pos, plotter,
                                                bindung_verschiebung)
        if endpunkt_liste is None:
            endpunkt_liste = []

        for endpunkt in endpunkt_liste:
            koordinaten = endpunkt[0]
            position_kette = endpunkt[1]

            phenyl_mesh = pv.Polygon(
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


def alkohol_substituent(stamm_kette_punkte, alkohol_input, plotter, bindung_verschiebung, verschiebung_h=0.2):
    try:
        alle_alkohol_alle_pos = dic_converter(alkohol_input)

        if alle_alkohol_alle_pos is None:
            alle_alkohol_alle_pos = {}

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_alkohol_alle_pos, plotter,
                                                bindung_verschiebung)
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

            wasserstoff_verbindung_linie = pv.lines_from_points([wasserstoff_anfangspunkt, wasserstoff_endpunkt])
            plotter.add_mesh(wasserstoff_verbindung_linie, line_width=3)
            besetzt_liste.append(position_kette)

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


def amino_substituent(stamm_kette_punkte, amin_input, plotter,
                      bindung_verschiebung, verschiebung_h=0.2):
    try:
        alle_amin_alle_pos = dic_converter(amin_input)

        if alle_amin_alle_pos is None:
            alle_amin_alle_pos = {}

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_amin_alle_pos, plotter, bindung_verschiebung)

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
            besetzt_liste.append(position_kette)

            vorzeichen = -1 if position_kette in besetzt_liste else 1
            y_formel = bindung_verschiebung - 2 if position_kette % 2 == 0 else bindung_verschiebung

            wasserstoff_verbindung_punkte_links = np.array([
                koordinaten,
                np.array([stamm_kette_punkte[position_kette][0] + vorzeichen * -bindung_verschiebung,
                          y_formel - verschiebung_h if position_kette % 2 == 0 else y_formel + verschiebung_h, 0])
            ])

            wasserstoff_verbindung_punkte_rechts = np.array([
                koordinaten,
                np.array([2 * koordinaten[0] - wasserstoff_verbindung_punkte_links[1][0],
                          y_formel - verschiebung_h if position_kette % 2 == 0 else y_formel + verschiebung_h, 0])
            ])

            wasserstoff_verbindung_linie_links = pv.lines_from_points(wasserstoff_verbindung_punkte_links)
            wasserstoff_verbindung_linie_rechts = pv.lines_from_points(wasserstoff_verbindung_punkte_rechts)

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


def saeure_substituent(stamm_kette_punkte, saeure_input, plotter, bindung_verschiebung, verschiebung_h=0.2):
    try:
        alle_saeure_alle_pos = dic_converter(saeure_input)
        zaehl_wort = saeure_input[0][1]
        if zaehl_wort == "di":
            alle_saeure_alle_pos["säure"] = [1, len(stamm_kette_punkte)]
        else:
            alle_saeure_alle_pos["säure"] = [1]

        endpunkt_liste = substituent_verbindung(stamm_kette_punkte, alle_saeure_alle_pos, plotter, bindung_verschiebung)
        verbindung_koordinaten = [endpunkt_liste[0][0], stamm_kette_punkte[0]]
        print(verbindung_koordinaten)
        bindung_zeichnen(verbindung_koordinaten, plotter, alle_bindungen_alle_pos={"en": (1, 2)},
                         verschiebung_bindung=-0.05,
                         laenge_bindung=0.008)


    except Exception as e:
        print(f"Fehler in der Darstellung der Säure / des Aldehyds: {e}")


def darsteller(stereo, alkan_input, halogen_input, phenyl_input, alkohol_input, amin_input, is_cyclo, stamm_input,
               bindung_typ,
               saeure_input, aldehyd_input, plotter, bindung_verschiebung):
    try:
        plotter.clear()
        bindung_verschiebung = -np.deg2rad(bindung_verschiebung) - (np.pi / 2)
        stamm_kette_punkte = stamm_kette(stamm_input)

        if stamm_kette_punkte is not None:
            alle_bindungen_alle_pos = dic_converter(bindung_typ)
            if alle_bindungen_alle_pos is not None:
                stamm_anpassung(stamm_kette_punkte, plotter, alle_bindungen_alle_pos)
                bindung_zeichnen(stamm_kette_punkte, plotter, alle_bindungen_alle_pos)

        else:
            print("Fehler in der Stammkette, setze fort")

        substituenten = [
            (alkan_substituent, alkan_input),
            (halogen_substituent, halogen_input),
            (phenyl_substituent, phenyl_input),
            (alkohol_substituent, alkohol_input),
            (amino_substituent, amin_input),
            (saeure_substituent, saeure_input)
        ]

        for funktion, text in substituenten:
            if text:
                funktion(stamm_kette_punkte, text, plotter, bindung_verschiebung)

        for kohlenstoff in stamm_kette_punkte:
            punkt = pv.Sphere(radius=0.04, center=kohlenstoff)
            plotter.add_mesh(punkt, line_width=4, color=(122, 20, 122))

        stamm_kette_zeichnung = pv.lines_from_points(stamm_kette_punkte)
        plotter.add_mesh(stamm_kette_zeichnung, line_width=4, color=(0, 0, 0))
        plotter.add_axes()
        plotter.camera_position = "xy"
        plotter.render()

    except Exception as e:
        print(f"Fehler im Darsteller: {e}")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Molekül Simulation")
        self.setGeometry(100, 100, 2200, 2000)

        layout_vertikal = QVBoxLayout()
        layout_horizontal = QHBoxLayout()

        self.textbox = QTextEdit(self)
        layout_horizontal.addWidget(self.textbox, stretch=3)
        self.textbox.setFont(QFont('Gill Sans MT', 15))

        self.dial = QDial()
        self.dial.setMinimum(0)
        self.dial.setMaximum(360)
        self.dial.setNotchesVisible(True)
        self.dial.setWrapping(True)
        layout_horizontal.addWidget(self.dial, stretch=1)

        layout_horizontal.setSpacing(10)

        layout_vertikal.addLayout(layout_horizontal)

        self.plotter = pvqt.QtInteractor(self)
        layout_vertikal.addWidget(self.plotter.interactor)

        self.button = QPushButton("Start")
        layout_vertikal.addWidget(self.button)
        self.button.clicked.connect(self.update_plot)

        self.dial.valueChanged.connect(self.update_plot)

        central_widget = QWidget()
        central_widget.setLayout(layout_vertikal)
        self.setCentralWidget(central_widget)

    def update_plot(self):
        input_text = self.textbox.toPlainText()
        print(input_text)
        resultat = text_auslesen(input_text)

        if not resultat:
            print("Keine Eingabe")
            return

        (stereo, alkan_input, halogen_input, phenyl_input, alkohol_input, amin_input, is_cyclo, stamm, bindung_typ,
         saeure_input, aldehyd_input) = resultat
        print(resultat)

        darsteller(
            stereo,
            alkan_input,
            halogen_input,
            phenyl_input,
            alkohol_input,
            amin_input,
            is_cyclo,
            stamm,
            bindung_typ,
            saeure_input,
            aldehyd_input,
            self.plotter,
            self.dial.value()
        )


if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()
