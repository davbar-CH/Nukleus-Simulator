import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtWidgets import *
import numpy as np
import time
import re
from Elemente import periodensystem
import random

def text_auslesen(input_text):
    """
    Die Funktion liest mit Regex den Text aus und berechnet die Neutronenzahl

    :param input_text: Text aus der Textbox, in der Form "Element-Massenzahl".
    :return: gibt den Namen des Elementes, die Massenzahl, die Neutronenzahl und Ordnungszahl zurück.
    """
    element = re.search(r"(\b[A-Z][a-z]*\b)", input_text)
    massenzahl = re.search(r"(\b[0-9]+\b)", input_text)

    if not element or not massenzahl:
        return None

    element = element.group()
    massenzahl = int(massenzahl.group())
    name, ordnungszahl, n_schalen = periodensystem[element]
    n_neutronen = massenzahl - int(ordnungszahl)

    print(f"{name}, Massenzahl:{massenzahl}, Ordnungszahl:{ordnungszahl}")
    return name, massenzahl, n_neutronen, ordnungszahl, n_schalen


def positions_berechnung(massenzahl, n_neutronen, ordnungszahl, n_schalen, plotter, radius=1):
    """
    Die Funktion berechnet anhand der Massenzahl (also der Anzahl Nukleonen) die Position der Nukleonen
    mit Hilfe des "Fibonacci Sphere Algorithm".

    Danach nähert die Funktion die Nukleonen (also Sphären) so nahe aneinander, bis sie sich berühren. In Blau werden
    die Neutronen, in Rot die Protonen dargestellt.

    Die Funktion berechnet zudem die Lage und Ausrichtung der Kreisbahnen der Elektronen (k, l, m usw. Schalen).

    :param radius: der Radius der einzelnen Nukleonen
    :param massenzahl: die Massenzahl des Elementes
    :param n_neutronen: die Anzahl Neutronen eines Elementes, wird in text_auslesen() berechnet
    :param ordnungszahl: die Ordnungszahl des Elementes
    :param plotter: das plotter-widget, um den Nukleus in pyvista darzustellen
    :return: startet automatisch den Plot
    """

    indices = np.arange(0, massenzahl, dtype=float) + 0.5
    indices_bahnen = np.arange(0, n_schalen, dtype=float) + 0.5

    phi = np.arccos(1 - 2 * indices / massenzahl)
    theta = np.pi * (1 + 5 ** 0.5) * indices

    phi_schalen = np.arccos(1 - 2 * indices_bahnen / n_schalen)
    theta_schalen = np.pi * (1 + 5 ** 0.5) * indices_bahnen

    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    x_schalen = np.cos(theta_schalen) * np.sin(phi_schalen)
    y_schalen = np.sin(theta_schalen) * np.sin(phi_schalen)
    z_schalen = np.cos(phi_schalen)


    start2 = time.time()
    begrenzung = radius * np.sqrt(massenzahl / np.pi) # oder 2r=2Rsin(n/π​)

    plotter.clear()

    alle_kugeln = []
    groesste_x_cord_2 = 0
    for x_cords, y_cords, z_cords in zip(x, y, z):
        p = np.array([x_cords, y_cords, z_cords])
        norm_p = np.linalg.norm(p)

        k = 1 - (begrenzung / norm_p)
        op_final = (1 - k) * p

        groesste_x_cord_1 = op_final[0]
        if groesste_x_cord_1 > groesste_x_cord_2:
            groesste_x_cord_2 = groesste_x_cord_1

        sphere = pv.Sphere(center=op_final, radius=radius).triangulate()
        alle_kugeln.append(sphere)

    abstand = 7

    for x_cords_schalen, y_cords_schalen, z_cords_schalen in zip(x_schalen, y_schalen, z_schalen):

        bahn = pv.Disc(c_res=50, inner=groesste_x_cord_2 + abstand, outer=groesste_x_cord_2 + abstand+0.5,
                       normal=(x_cords_schalen, y_cords_schalen, z_cords_schalen))
        plotter.add_mesh(bahn)
        abstand = abstand + 3


    if massenzahl:
        random_neutronen = random.sample(alle_kugeln, n_neutronen)

        for element in random_neutronen:
            alle_kugeln.remove(element)

        random_protonen = random.sample(alle_kugeln, ordnungszahl)

        for neutron in random_neutronen:
            plotter.add_mesh(neutron, color="blue", style='wireframe')

        for proton in random_protonen:
            plotter.add_mesh(proton, color="red", style='wireframe')

        plotter.add_axes()
        plotter.render()

    end2 = time.time()
    print(f"Laufzeit: {end2 - start2}")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Nukleus Simulation")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout()
        central_widget.setLayout(layout)

        self.textbox = QTextEdit(self)
        layout.addWidget(self.textbox)

        self.plotter = pvqt.QtInteractor(self)
        layout.addWidget(self.plotter.interactor)

        self.button = QPushButton("Start")
        layout.addWidget(self.button)
        self.button.clicked.connect(self.update_plot)


    def update_plot(self):
        input_text = self.textbox.toPlainText()

        resultat = text_auslesen(input_text)

        if not resultat:
            print("Keine Eingabe")
            return

        name, massenzahl, n_neutronen, ordnungszahl, n_schalen = resultat

        positions_berechnung(
            massenzahl,
            n_neutronen,
            ordnungszahl,
            n_schalen,
            self.plotter
        )


app = QApplication([])
window = MainWindow()
window.show()
app.exec_()
