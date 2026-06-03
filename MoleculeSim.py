import pyvista as pv
import pyvistaqt as pvqt
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
    |phenyl|hydroxy|amino|oxo|formyl)"

    stereo = re.findall(stereo_pattern, input_text, flags=re.IGNORECASE)

    halogen = re.findall(pattern_halogen, input_text, flags=re.IGNORECASE)
    alkan = re.findall(pattern_alkan, input_text, flags=re.IGNORECASE)
    phenyl = re.findall(pattern_phenyl, input_text, flags=re.IGNORECASE)
    amin = re.findall(amin_pattern, input_text, flags=re.IGNORECASE)
    alkohol = re.findall(alkohol_pattern, input_text, flags=re.IGNORECASE)

    input_ohne_stereo = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)
    input_ohne_stereo_sub = re.sub(alle_sub_pattern, '', input_ohne_stereo, flags=re.IGNORECASE)

    isCyclo = False
    if re.search(r"cyclo", input_ohne_stereo_sub, re.IGNORECASE):
        isCyclo = True

    stamm_pattern = r"(?:cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
    stamm = re.findall(stamm_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    bindung_pattern = rf"{zaehl_woerter_pattern}(en|in)(?!\w)"
    bindung_typ = re.findall(bindung_pattern, input_ohne_stereo_sub, flags=re.IGNORECASE)

    endung_pattern_saeure_al = rf"{zaehl_woerter_pattern}(säure|al)"
    endung_saeure_al = re.findall(endung_pattern_saeure_al, input_ohne_stereo_sub, flags=re.IGNORECASE)

    return stereo, alkan, halogen, phenyl, isCyclo, stamm, bindung_typ, endung_saeure_al, amin, alkohol


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


def stamm_kette(stamm, plotter):
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

        stamm = stamm[0].lower()
        stamm_kette_punkte = np.array([[x * 0.5, (1 - (-1) ** x) / 2, 0] for x in range(0, stamm_laenge.get(stamm))])

        for kohlenstoff in stamm_kette_punkte:
            punkt = pv.Sphere(radius=0.04, center=kohlenstoff)
            plotter.add_mesh(punkt, line_width=4, color=(122, 20, 122))

        stamm_kette = pv.lines_from_points(stamm_kette_punkte)
        plotter.add_mesh(stamm_kette, line_width=4, color=(0, 0, 0))

        return stamm_kette_punkte

    except Exception as e:
        print(f"Fehler in der stamm_kette: {e}")


def stamm_anpassung(stamm_kette_punkte, plotter, bindung_typ):
    try:
        alle_bindungen_alle_pos = dic_converter(bindung_typ)

        alle_positionen_en = alle_bindungen_alle_pos.get("en")
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

        alle_positionen_in = alle_bindungen_alle_pos.get("in")
        if alle_bindungen_alle_pos:
            for pos in alle_positionen_in:

                if pos - 1 == 0:
                    new_val = 0

                else:
                    new_val = 1 - stamm_kette_punkte[pos - 1][1]

                    for j in range(pos + 1, len(stamm_kette_punkte)):
                        stamm_kette_punkte[j][1] = new_val if j == pos + 1 else 1 - stamm_kette_punkte[j][1]

                stamm_kette_punkte[pos - 1][1] = new_val
                stamm_kette_punkte[pos][1] = new_val

        else:
            print("Keine Bindung")

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


def alkan_substituent(stamm_kette_punkte, alkan_input, plotter):
    try:
        """
        pos = Position des Substituents, in der Form 1 oder 3,8
        sub = Substituent, in der Form Methyl oder methyl

        Ist eine Position besetzt, dann wird das Vorzeichen gekehrt, damit die Substituenten trotzdem angezeigt werden können

        Bei ungeraden Positionen wird bei y=0 gestartet, bei geraden Positionen bei y=1

        Der Substituent startet bei der jeweiligen x-Koordinate (pos-1). Die x-Koordinate alterniert jeweils zwischen 
        der x-Koordinate und x-Koordinate - 0.25. Die y-Koordinate wird stets um 0.5 grösser (bei geraden Positionen) bzw.
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

        global besetzt_liste
        besetzt_liste = []

        for alkan in alle_alkane_alle_pos:
            for sub_pos in alle_alkane_alle_pos.get(alkan):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = (lambda x: (x + 2) * 0.5) if sub_pos % 2 == 0 else (lambda x: x * -0.5)
                laenge = substituent_laenge.get(alkan) + 1

                sub_alkan_punkte = np.array([
                    [stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * ((-1 + (-1) ** x) / 8), y_formel(x), 0]
                    for x in range(0, laenge)
                ])

                besetzt_liste.append(sub_pos)
                alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                plotter.add_mesh(alkan_kette, line_width=2)

    except:
        print("kein Substituent")


def halogen_substituent(stamm_kette_punkte, halogen_input, plotter):
    try:
        # Name, Farbe, Grösse
        halogen_zeichnung = {
            "fluor": ('F', "#FFD1DC", 40),
            "chlor": ('Cl', "#228B22", 50),
            "brom": ('Br', "#CC5500", 60),
            "iod": ('I', "#9D00FF", 70)
        }

        alle_halogene_alle_pos = dic_converter(halogen_input)

        for halogen in alle_halogene_alle_pos:
            for sub_pos in alle_halogene_alle_pos.get(halogen):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

                halogen_verbindung_punkte = np.array([
                    np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
                ])

                besetzt_liste.append(sub_pos)
                verbindung_halogen = pv.lines_from_points(halogen_verbindung_punkte)
                plotter.add_mesh(verbindung_halogen, line_width=2)

                halogen_lower = halogen.lower()
                zeichnung_pos = np.array([[stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25,y_formel, 0]])

                text, color, point_size = halogen_zeichnung[halogen_lower]

                plotter.add_point_labels(
                    points=zeichnung_pos,
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


def phenyl_substituent(stamm_kette_punkte, phenyl_input, plotter):
    try:
        # phenyl | hydroxy | amino | oxo | formyl
        alle_phenyl_sub_alle_pos = dic_converter(phenyl_input)

        for phenyl in alle_phenyl_sub_alle_pos:
            for sub_pos in alle_phenyl_sub_alle_pos.get(phenyl):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

                # die nachfolgenden Zahlen sind zufällig gewählt, sieht einfach am besten aus
                phenyl_verbindung_punkte = np.array([
                    np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
                ])
                endpunkt = phenyl_verbindung_punkte[1]

                besetzt_liste.append(sub_pos)
                verbindung_phenyl = pv.lines_from_points(phenyl_verbindung_punkte)
                plotter.add_mesh(verbindung_phenyl, line_width=2)

                phenyl_mesh = pv.Polygon(
                    center=(endpunkt[0], endpunkt[1] + 0.5 if sub_pos % 2 == 0 else endpunkt[1] - 0.5, 0),
                    radius=0.5,
                    fill=False)
                plotter.add_mesh(phenyl_mesh, line_width=2, color=(0, 0, 0))
                punkte = phenyl_mesh.points

                bindung_zeichnen(punkte, plotter, alle_bindungen_alle_pos={"en": (1, 3, 5)}, verschiebung_bindung=-0.05,
                        laenge_bindung=0.008)

    except Exception as e:
        print(f"Fehler in der Darstellung von Phenyl: {e}")

def alkohol_substituent(stamm_kette_punkte, alkohol_input, plotter, verschiebung_H=0.2):
    try:
        alle_alkohol_alle_pos = dic_converter(alkohol_input)

        for alkohol in alle_alkohol_alle_pos:
            for sub_pos in alle_alkohol_alle_pos.get(alkohol):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

                # die nachfolgenden Zahlen sind zufällig gewählt, sieht einfach am besten aus
                alkohol_verbindung_punkte = np.array([
                    np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
                ])
                endpunkt_verbindung = alkohol_verbindung_punkte[1]

                plotter.add_point_labels(
                    points=endpunkt_verbindung,
                    labels=["O"],
                    font_size=40,
                    point_color="#ec0c0d",
                    point_size=40,
                    render_points_as_spheres=True,
                    always_visible=True,
                    shape=None
                )
                besetzt_liste.append(sub_pos)
                verbindung_alkohol = pv.lines_from_points(alkohol_verbindung_punkte)
                plotter.add_mesh(verbindung_alkohol, line_width=2)

                wasserstoff_verbindung_punkte = np.array([
                    endpunkt_verbindung,
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25,
                              y_formel + verschiebung_H if sub_pos % 2 == 0 else y_formel - verschiebung_H, 0])
                ])

                wasserstoff_verbindung_linie = pv.lines_from_points(wasserstoff_verbindung_punkte)
                plotter.add_mesh(wasserstoff_verbindung_linie, line_width=3)

                plotter.add_point_labels(
                    points=[wasserstoff_verbindung_punkte[1]],
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

def amino_substituent(stamm_kette_punkte, amin_input, plotter, verschiebung_H=0.2):
    try:
        alle_amin_alle_pos = dic_converter(amin_input)

        for amin in alle_amin_alle_pos:
            for sub_pos in alle_amin_alle_pos.get(amin):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

                amin_verbindung_punkte = np.array([
                    np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
                ])
                endpunkt_verbindung = amin_verbindung_punkte[1]

                plotter.add_point_labels(
                    points=endpunkt_verbindung,
                    labels=["N"],
                    font_size=40,
                    point_color="#1e7fcb",
                    point_size=40,
                    render_points_as_spheres=True,
                    always_visible=True,
                    shape=None
                )
                besetzt_liste.append(sub_pos)
                verbindung_amin = pv.lines_from_points(amin_verbindung_punkte)
                plotter.add_mesh(verbindung_amin, line_width=2)

                wasserstoff_verbindung_punkte_links = np.array([
                    endpunkt_verbindung,
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25,
                              y_formel + verschiebung_H if sub_pos % 2 == 0 else y_formel - verschiebung_H, 0])
                ])

                wasserstoff_verbindung_punkte_rechts = np.array([
                    endpunkt_verbindung,
                    np.array([2 * endpunkt_verbindung[0] - wasserstoff_verbindung_punkte_links[1][0],
                         y_formel + verschiebung_H if sub_pos % 2 == 0 else y_formel - verschiebung_H, 0])
                ])

                wasserstoff_verbindung_linie_links = pv.lines_from_points(wasserstoff_verbindung_punkte_links)
                wasserstoff_verbindung_linie_rechts = pv.lines_from_points(wasserstoff_verbindung_punkte_rechts)

                plotter.add_mesh(wasserstoff_verbindung_linie_links, line_width=3)
                plotter.add_mesh(wasserstoff_verbindung_linie_rechts, line_width=3)

                plotter.add_point_labels(
                    points=[wasserstoff_verbindung_punkte_links[1], wasserstoff_verbindung_punkte_rechts[1]],
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

def darsteller(stereo, alkan_input, halogen_input, phenyl_input, isCyclo, stamm, bindung_typ,
               endung_saeure_al, alkohol_input, amin_input, plotter):
    try:
        plotter.clear()

        stamm_kette_punkte = stamm_kette(stamm, plotter)
        alle_bindungen_alle_pos = stamm_anpassung(stamm_kette_punkte, plotter, bindung_typ)
        bindung_zeichnen(stamm_kette_punkte, plotter, alle_bindungen_alle_pos)

        substituenten = [
            (alkan_substituent, alkan_input),
            (halogen_substituent, halogen_input),
            (phenyl_substituent, phenyl_input),
            (alkohol_substituent, alkohol_input),
            (amino_substituent, amin_input),
        ]

        for funktion, text in substituenten:
            if text:
                funktion(stamm_kette_punkte, text, plotter)

        plotter.add_axes()
        plotter.camera_position = "xy"
        plotter.render()

    except Exception as e:
        print(f"Fehler im Darsteller: {e}")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Molekül Simulation")
        self.setGeometry(100, 100, 800, 600)

        layout_vertikal = QVBoxLayout()
        layout_horizontal = QHBoxLayout()

        self.infos = QLabel("")
        layout_horizontal.addWidget(self.infos)

        self.textbox = QTextEdit(self)
        layout_horizontal.addWidget(self.textbox)
        self.textbox.setFontPointSize(12)

        layout_horizontal.setSpacing(150)

        layout_vertikal.addLayout(layout_horizontal)

        self.plotter = pvqt.QtInteractor(self)
        layout_vertikal.addWidget(self.plotter.interactor)

        self.button = QPushButton("Start")
        layout_vertikal.addWidget(self.button)
        self.button.clicked.connect(self.update_plot)

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

        stereo, alkan_input, halogen_input, phenyl_input, isCyclo, stamm, bindung_typ, endung_saeure_al, amin_input, alkohol_input = resultat
        print(resultat)

        darsteller(
            stereo,
            alkan_input,
            halogen_input,
            phenyl_input,
            isCyclo,
            stamm,
            bindung_typ,
            endung_saeure_al,
            alkohol_input,
            amin_input,
            self.plotter
        )


if __name__ == "__main__":
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()
