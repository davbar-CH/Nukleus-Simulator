from pyvista import lines_from_points

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

        stamm_kette = lines_from_points(stamm_kette_punkte)
        plotter.add_mesh(stamm_kette, line_width=4, color=(0, 0, 0))

        return alle_bindungen_alle_pos

    except Exception as e:
        print(f"Fehler in der stamm_anpassung: {e}")