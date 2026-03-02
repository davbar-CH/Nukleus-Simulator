import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtWidgets import *
import numpy as np
import time
import re
from Elemente import periodensystem
import random

def text_auslesen(input_text):
    element = re.search(r"(\b[A-Z][a-z]*\b)", input_text)
    massenzahl = re.search(r"(\b[0-9]+\b)", input_text)

    if not element or not massenzahl:
        return None

    element = element.group()
    massenzahl = int(massenzahl.group())
    name, ordnungszahl = periodensystem[element]
    n_neutronen = massenzahl - int(ordnungszahl)

    print(name, massenzahl, n_neutronen)
    return name, massenzahl, n_neutronen, ordnungszahl


def positions_berechnung(massenzahl, n_neutronen, ordnungszahl, plotter):
    radius = 1
    indices = np.arange(0, massenzahl, dtype=float) + 0.5

    phi = np.arccos(1 - 2 * indices / massenzahl)
    theta = np.pi * (1 + 5 ** 0.5) * indices

    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)

    start2 = time.time()
    begrenzung = radius * np.sqrt(massenzahl / np.pi)


    alle_kugeln = []
    for x_cords, y_cords, z_cords in zip(x, y, z):
        p = np.array([x_cords, y_cords, z_cords])
        norm_p = np.linalg.norm(p)

        k = 1 - (begrenzung / norm_p)
        op_final = (1 - k) * p

        sphere = pv.Sphere(center=op_final, radius=radius).triangulate()
        alle_kugeln.append(sphere)

    plotter.clear()

    if massenzahl:
        random_neutronen = random.sample(alle_kugeln, n_neutronen)

        for element in random_neutronen:
            alle_kugeln.remove(element)

        random_protonen = random.sample(alle_kugeln, ordnungszahl)

        for neutron in random_neutronen:
            plotter.add_mesh(neutron, color="blue", style='wireframe')

        for proton in random_protonen:
            plotter.add_mesh(proton, color="red", style='wireframe')

        plotter.render()

    end2 = time.time()
    print(end2 - start2)


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

        name, massenzahl, n_neutronen, ordnungszahl = resultat

        positions_berechnung(
            massenzahl,
            n_neutronen,
            ordnungszahl,
            self.plotter
        )


app = QApplication([])
window = MainWindow()
window.show()
app.exec_()