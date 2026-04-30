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
    substituenten_pattern = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                             r"fluor|chlor|brom|iod|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|phenyl|nitro"
                             r"|hydroxy|amino|oxo)")

    stereo = re.findall(stereo_pattern, input_text, flags=re.IGNORECASE)

    substituent = re.findall(substituenten_pattern, input_text, flags=re.IGNORECASE)

    it_ohne_ster = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)

    it_ohne_ster_sub = re.sub(substituenten_pattern, '', it_ohne_ster, flags=re.IGNORECASE)

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

    return stereo, substituent, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest


def stamm_kette(stamm, plotter):
    try:
        stamm_laenge = {
            "Eth": 2,
            "eth": 2,
            "Prop": 3,
            "prop": 3,
            "But": 4,
            "but": 4,
            "Pent": 5,
            "pent": 5,
            "Hex": 6,
            "hex": 6,
            "Hept": 7,
            "hept": 7,
            "Oct": 8,
            "oct": 8,
            "Non": 9,
            "non": 9,
            "Dec": 10,
            "dec": 10,
        }

        stamm_kette_punkte = np.array([[x * 0.5, (1 - (-1) ** x) / 2, 0] for x in range(0, stamm_laenge.get(stamm[0]))])
        stamm_kette = pv.lines_from_points(stamm_kette_punkte)
        plotter.add_mesh(stamm_kette, line_width=4,color=(0,0,0))
        return stamm_kette_punkte

    except Exception as e:
        print(f"Fehler in der stamm_kette: {e}")

def bindung(stamm_kette_punkte, bindung_typ, plotter):
    alle_pos = [x[0] for x in bindung_typ]
    alle_bindungen = [x[2] for x in bindung_typ]
    alle_bindungen_alle_pos = {}

    verschiebung_oben = np.array([0,0.3,0])
    verschiebung_unten_links = np.array([-0.24, -0.18,0])
    verschiebung_unten_rechts = np.array([0.24, -0.18, 0])

    for i, position in enumerate(alle_pos):
        bindung_pos = {alle_bindungen[i]: [int(x) for x in re.findall(r"\d", position)]}
        alle_bindungen_alle_pos.update(bindung_pos)

    for bindung in alle_bindungen_alle_pos:
        for bindung_pos in alle_bindungen_alle_pos.get(bindung):

            if bindung_pos % 2 == 1:
                sub_alken_punkte = np.array([
                    verschiebung_oben+np.array([stamm_kette_punkte[bindung_pos - 1][0], stamm_kette_punkte[bindung_pos - 1][1], 0]),
                    verschiebung_unten_links+np.array([stamm_kette_punkte[bindung_pos][0], stamm_kette_punkte[bindung_pos][1], 0])
                ])
                if bindung == "in":
                    sub_alkin_punkte = np.array([
                        verschiebung_oben + np.array(
                            [stamm_kette_punkte[bindung_pos - 1][0], stamm_kette_punkte[bindung_pos - 1][1], 0]),
                        verschiebung_unten_links + np.array(
                            [stamm_kette_punkte[bindung_pos][0], stamm_kette_punkte[bindung_pos][1], 0])
                    ])


            else:
                sub_alken_punkte = np.array([
                    verschiebung_unten_rechts + np.array(
                        [stamm_kette_punkte[bindung_pos - 1][0], stamm_kette_punkte[bindung_pos - 1][1], 0]),
                    verschiebung_oben + np.array(
                        [stamm_kette_punkte[bindung_pos][0], stamm_kette_punkte[bindung_pos][1], 0])
                ])

            sub_alkan_kette = pv.lines_from_points(sub_alken_punkte)
            plotter.add_mesh(sub_alkan_kette, line_width=2, color=(255, 0, 0))


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
            "Methyl": 1,
            "methyl": 1,
            "Ethyl": 2,
            "ethyl": 2,
            "Propyl": 3,
            "propyl": 3,
            "Butyl": 4,
            "butyl": 4,
            "Pentyl": 5,
            "pentyl": 5,
            "Hexyl": 6,
            "hexyl": 6,
            "Heptyl": 7,
            "heptyl": 7,
            "Octyl": 8,
            "octyl": 8,
            "Nonyl": 9,
            "nonyl": 9,
            "Decyl": 10,
            "decyl": 10,
        }

        alle_pos = [x[0] for x in substituent]
        alle_sub = [x[2] for x in substituent]
        alle_sub_alle_pos = {}

        for i, position in enumerate(alle_pos):
            sub_pos = {alle_sub[i]: [int(x) for x in re.findall(r"\d", position)]}
            alle_sub_alle_pos.update(sub_pos)

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


def darsteller(stereo, substituent, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest, plotter):
    plotter.clear()

    stamm_kette_punkte = stamm_kette(stamm, plotter)
    bindung(stamm_kette_punkte, bindung_typ, plotter)
    alkan_substituent(stamm_kette_punkte, substituent, plotter)

    plotter.add_axes()
    plotter.camera_position = "xy"
    plotter.render()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Nukleus Simulation")
        self.setGeometry(100, 100, 800, 600)

        layout_vertikal = QVBoxLayout()
        layout_horizontal = QHBoxLayout()

        self.infos = QLabel("")
        layout_horizontal.addWidget(self.infos)

        self.textbox = QTextEdit(self)
        layout_horizontal.addWidget(self.textbox)

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

        stereo, substituent, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest = resultat
        print(stereo, substituent, isCyclo, stamm, bindung_typ, endung_saeure_al, endung_rest)

        darsteller(
            stereo,
            substituent,
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
