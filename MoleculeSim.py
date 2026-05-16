import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtWidgets import *
import numpy as np
import re


def text_auslesen(input_text):
    """
    Die Funktion liest mit Regex den Text aus und berechnet die Neutronenzahl

    :param input_text: Text aus der Textbox, in der Form "Element-Massenzahl".
    :return: gibt den Namen des Elementes, die Massenzahl, die Neutronenzahl und Ordnungszahl zurück.
    """
    stereo_pattern = r"\(([EZRSezrs,\d]+)\)-"
    substituenten_pattern_halogen = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                                     r"fluor|chlor|brom|iod)")

    substituenten_pattern_alkan = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                                   r"methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|octyl|nonyl|decyl)")

    substituenten_pattern_rest = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                                  r"phenyl|hydroxy|amino|oxo|formyl)")

    alle_sub_pattern = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                        r"fluor|chlor|brom|iod|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|phenyl"
                        r"|hydroxy|amino|oxo|formyl)")

    stereo = re.findall(stereo_pattern, input_text, flags=re.IGNORECASE)

    substituent_alkan = re.findall(substituenten_pattern_alkan, input_text, flags=re.IGNORECASE)

    substituent_halogen = re.findall(substituenten_pattern_halogen, input_text, flags=re.IGNORECASE)

    substituent_rest = re.findall(substituenten_pattern_rest, input_text, flags=re.IGNORECASE)

    it_ohne_ster = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)

    it_ohne_ster_sub = re.sub(alle_sub_pattern, '', it_ohne_ster, flags=re.IGNORECASE)

    isCyclo = False
    if re.search(r"cyclo", it_ohne_ster_sub, re.IGNORECASE):
        isCyclo = True

    stamm_pattern = r"(?:cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
    stamm = re.findall(stamm_pattern, it_ohne_ster_sub, flags=re.IGNORECASE)

    bindung_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?(en|in)(?!\w)"
    bindung_typ = re.findall(bindung_pattern, it_ohne_ster_sub, flags=re.IGNORECASE)

    endung_pattern_saeure_al = r"(säure|al)"
    endung_saeure_al = re.findall(endung_pattern_saeure_al, it_ohne_ster_sub, flags=re.IGNORECASE)

    endung_pattern_rest = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?(ol|amin)"
    endung_rest = re.findall(endung_pattern_rest, it_ohne_ster_sub, flags=re.IGNORECASE)

    return stereo, substituent_alkan, substituent_halogen, substituent_rest, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest


def dic_converter(input):
    alle_pos = [x[0] for x in input]
    alle_gruppen = [x[2].lower() for x in input]

    alle_gruppen_alle_pos = {}
    for i, position in enumerate(alle_pos):
        bindung_pos = {alle_gruppen[i]: [int(x) for x in re.findall(r"\d", position)]}
        alle_gruppen_alle_pos.update(bindung_pos)

    return alle_gruppen_alle_pos


def stamm_kette(stamm, plotter, bindung_typ):
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

        alle_bindungen_alle_pos = dic_converter(bindung_typ)
        alle_positionen_en = alle_bindungen_alle_pos.get("en")

        if "en" in alle_bindungen_alle_pos:
            for i, pos in enumerate(alle_positionen_en):

                if i + 1 == len(alle_positionen_en):
                    break

                else:
                    if alle_positionen_en[i + 1] - pos == 1:
                        if stamm_kette_punkte[pos][1] == 0:
                            stamm_kette_punkte[pos][1] = 1
                        elif stamm_kette_punkte[pos][1] == 1:
                            stamm_kette_punkte[pos][1] = 0
                    else:
                        pass

        alle_positionen_in = alle_bindungen_alle_pos.get("in")
        if "in" in alle_bindungen_alle_pos:
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

        for kohlenstoff in stamm_kette_punkte:
            punkt = pv.Sphere(radius=0.04, center=kohlenstoff)
            plotter.add_mesh(punkt, line_width=4, color=(122, 20, 122))

        return stamm_kette_punkte, alle_bindungen_alle_pos

    except Exception as e:
        print(f"Fehler in der stamm_kette: {e}")


def bindung(stamm_kette_punkte, plotter, alle_bindungen_alle_pos, verschiebung_bindung=0.1, laenge_bindung=0.2):
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


def alkan_substituent(stamm_kette_punkte, substituent, plotter):
    try:
        """
        pos = Position des Substituents, in der Form 1 oder 3,8
        sub = Substituent, in der Form Methyl oder methyl

        ('2,3', 'Di', 'ethyl'), ('4,5', 'di', 'methyl') Form wird in {'ethyl': [2, 3], 'methyl': [5, 5]} umgewandelt
        Zahlen sind integer

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

        alle_sub_alle_pos = dic_converter(substituent)

        global besetzt_liste
        besetzt_liste = []

        for alkan_substituent in alle_sub_alle_pos:
            for sub_pos in alle_sub_alle_pos.get(alkan_substituent):
                vorzeichen = -1 if sub_pos in besetzt_liste else 1
                y_formel = (lambda x: (x + 2) * 0.5) if sub_pos % 2 == 0 else (lambda x: x * -0.5)
                laenge = substituent_laenge.get(alkan_substituent) + 1

                sub_alkan_punkte = np.array([
                    [stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * ((-1 + (-1) ** x) / 8), y_formel(x), 0]
                    for x in range(0, laenge)
                ])

                besetzt_liste.append(sub_pos)
                sub_alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                plotter.add_mesh(sub_alkan_kette, line_width=2)

    except:
        print("kein Substituent")


def element_substituent(stamm_kette_punkte, substituent_element, plotter):
    # Name, Farbe, Grösse
    element_zeichnung = {
        "wasserstoff":("H","#d9e4ea",20),
        "sauerstoff":("O","#ec0c0d", 40),
        "fluor": ('F', "#FFD1DC", 40),
        "chlor": ('Cl', "#228B22", 50),
        "brom": ('Br', "#CC5500", 60),
        "iod": ('I', "#9D00FF", 70)
    }

    alle_sub_alle_pos = dic_converter(substituent_element)

    for element in alle_sub_alle_pos:
        for sub_pos in alle_sub_alle_pos.get(element):
            vorzeichen = -1 if sub_pos in besetzt_liste else 1
            y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

            element_verbindung_punkte = np.array([
                np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
            ])

            besetzt_liste.append(sub_pos)
            verbindung_element = pv.lines_from_points(element_verbindung_punkte)
            plotter.add_mesh(verbindung_element, line_width=2)

            element_lower = element.lower()
            zeichnung_pos = np.array([[stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25,
                                       y_formel, 0]])

            text, color, point_size = element_zeichnung[element_lower]

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


def rest_substituent(stamm_kette_punkte, substituent_rest, plotter):
    #phenyl | hydroxy | amino | oxo | formyl
    alle_rest_sub_alle_pos = dic_converter(substituent_rest)

    for substituent in alle_rest_sub_alle_pos:
        for sub_pos in alle_rest_sub_alle_pos.get(substituent):
            vorzeichen = -1 if sub_pos in besetzt_liste else 1
            y_formel = 1.5 if sub_pos % 2 == 0 else -0.5

            # die nachfolgenden Zahlen sind zufällig gewählt, sieht einfach am besten aus
            if substituent == "phenyl":
                sub_verbindung_punkte = np.array([
                    np.array([stamm_kette_punkte[sub_pos - 1][0], stamm_kette_punkte[sub_pos - 1][1], 0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] + vorzeichen * -0.25, y_formel, 0])
                ])
                endpunkt = sub_verbindung_punkte[1]

                besetzt_liste.append(sub_pos)
                verbindung_substituent = pv.lines_from_points(sub_verbindung_punkte)
                plotter.add_mesh(verbindung_substituent, line_width=2)

                phenyl = pv.Polygon(center=(endpunkt[0], endpunkt[1] + 0.5 if sub_pos % 2 == 0 else endpunkt[1] - 0.5, 0), radius=0.5, fill=False)
                plotter.add_mesh(phenyl, line_width=2, color=(0, 0, 0))
                punkte = phenyl.points

                bindung(punkte, plotter, alle_bindungen_alle_pos={"en":(1,3,5)}, verschiebung_bindung=-0.05, laenge_bindung=0.008)

            if substituent == "hydroxy":
                substituent_element_sauerstoff = [(str(sub_pos), '', 'sauerstoff')] # der dic_converter nimmt die position nur als string entgegen
                element_substituent(stamm_kette_punkte, substituent_element_sauerstoff, plotter)

                wasserstoff_anfangspunkt = np.array([
                    np.array([0,0,0]),
                    np.array([stamm_kette_punkte[sub_pos - 1][0] - 0.25, y_formel, 0])
                ])

                substituent_element_wasserstoff = [(str(2), '', 'wasserstoff')]
                element_substituent(wasserstoff_anfangspunkt, substituent_element_wasserstoff, plotter)

            if substituent == "amino":
                pass

            if substituent == "oxo":
                pass

            if substituent == "formyl":
                pass

def darsteller(stereo, substituent_alkan, substituent_halogen, substituent_rest, isCyclo, stamm, bindung_typ,
               endung_saeure_al, endung_rest, plotter):
    plotter.clear()

    stamm_kette_punkte, alle_bindungen_alle_pos = stamm_kette(stamm, plotter, bindung_typ)
    bindung(stamm_kette_punkte, plotter, alle_bindungen_alle_pos)
    alkan_substituent(stamm_kette_punkte, substituent_alkan, plotter)
    element_substituent(stamm_kette_punkte, substituent_halogen, plotter)
    rest_substituent(stamm_kette_punkte, substituent_rest, plotter)

    plotter.add_axes()
    plotter.camera_position = "xy"
    plotter.render()


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

        stereo, substituent_alkan, substituent_halogen, substituent_rest, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest = resultat
        print(resultat)

        darsteller(
            stereo,
            substituent_alkan,
            substituent_halogen,
            substituent_rest,
            isCyclo,
            stamm,
            bindung_typ,
            endung_saeure_al,
            endung_rest,
            self.plotter
        )


app = QApplication([])
window = MainWindow()
window.show()
app.exec_()
