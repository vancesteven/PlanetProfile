import reaktoro as rkt
import os
import streamlit as st
#database_folder = "PlanetProfile/PlanetProfile/Thermodynamics/Reaktoro/Databases"

def read_salt_db(db_name):
    # Start from where SaltLoader.py is
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # Go up two directories, then into PlanetProfile/Thermodynamics/Reaktoro/Databases
    database_folder = os.path.abspath(
        os.path.join(base_dir, "..", "..", "PlanetProfile", "Thermodynamics", "Reaktoro", "Databases")
    )

    # If db_name has no extension, try adding ".dat"
    if not os.path.splitext(db_name)[1]:
        db_name += ".dat"


    db_path = os.path.join(database_folder, db_name)

    #st.success(f"Looking for databases at: {db_path}")
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database file does not exist at: {db_path}")

    aqueous_species = []
    try:
        db = rkt.PhreeqcDatabase.fromFile(db_path)
    except Exception:
        db = rkt.SupcrtDatabase.fromFile(db_path)

    for species in db.speciesWithAggregateState(rkt.AggregateState.Aqueous):
        aqueous_species.append(species.name())

    #print(aqueous_species)
    return aqueous_species
