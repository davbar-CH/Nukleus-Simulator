from pyvista import lines_from_points, Sphere
from pyvistaqt import QtInteractor
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import *
import numpy as np
from InputParser import *
from Substituenten import *
from ZweiDreiBindung import stamm_anpassung
from HelferGeometrie import bindung_zeichnen
from StammKette import stamm_kette

def darsteller(stereo, alkan_input, halogen_input, phenyl_input, alkohol_input, amin_input, is_cyclo, stamm_input,
               bindung_typ,
               saeure_input, aldehyd_input, plotter, bindung_verschiebung):
    try:
        plotter.clear()
        bindung_verschiebung = -np.deg2rad(bindung_verschiebung) - (np.pi / 2)
        stamm_kette_punkte, besetzt_liste = stamm_kette(stamm_input)

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
                funktion(stamm_kette_punkte, text, plotter, bindung_verschiebung, besetzt_liste)

        for kohlenstoff in stamm_kette_punkte:
            punkt = Sphere(radius=0.04, center=kohlenstoff)
            plotter.add_mesh(punkt, line_width=4, color=(122, 20, 122))

        stamm_kette_zeichnung = lines_from_points(stamm_kette_punkte)
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

        self.plotter = QtInteractor(self)
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

        (stereo, alkan_input, halogen_input, phenyl_input, alkohol_input, amin_input, is_cyclo, stamm_input, bindung_typ,
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
            stamm_input,
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
