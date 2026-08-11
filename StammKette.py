from numpy import array

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

        besetzt_liste = []
        stamm_kette_punkte = []

        stamm = stamm_input[0].lower()
        laenge = stamm_laenge.get(stamm)
        if laenge is not None:
            stamm_kette_punkte = array(
                [[x * dehnung_x, (1 - (-1) ** x) / dehnung_y, 0] for x in range(0, laenge)])
        else:
            print(f"Keine solche Stammkette{stamm}")

        return stamm_kette_punkte, besetzt_liste

    except Exception as e:
        print(f"Fehler in der stamm_kette: {e}")