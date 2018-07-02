#!/usr/bin/env python
#! -*- coding: utf-8 -*-

##    Description    eTOX repeat-dose toxicity extraction tool
##
##    Authors:       Elisabet Gregori (elisabet.gregori@upf.edu)
##                   Ignacio Pasamontes (ignacio.pasamontes@upf.edu)
##
##    Copyright 2018 Elisabet Gregori & Ignacio Pasamontes
##
##    studyExtraction is free software: you can redistribute it 
##    and/or modify it under the terms of the GNU General Public 
##    License as published by the Free Software Foundation version 3.
##
##    studyExtraction is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this code. If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import argparse, os
import psycopg2, pickle
from rdkit.Chem import PandasTools

def run(args):

    # conn = psycopg2.connect(host='172.20.16.76', dbname='vitic2016_1', 
    #                         user=args.user, password=args.passw)
    conn = psycopg2.connect(host='172.20.16.76', dbname='vitic2016_test', 
                            user=args.user, password=args.passw)
    curs = conn.cursor()

    #
    # Create and store compound dataframe
    #
    ## In the Postgres DB smiles are all lower case!!!! 
    # cmd = "SELECT smiles, subst_id, casno AS cas_number, common_name, molecular_weight, \
    #         pharmacological_action \
    #         FROM SUBSTANCE_IDS"    
    # df = pd.read_sql(cmd, con=conn)
    # df.pharmacological_action = df.pharmacological_action.str.capitalize()
    # # Add molecule object
    # smiles_df = pd.DataFrame(df.smiles.dropna().unique(), columns=['smiles'])
    # PandasTools.AddMoleculeColumnToFrame(smiles_df,'smiles','molecule',
    #                                     includeFingerprints=True)
    # df = pd.merge(df, smiles_df, 'left', on='smiles')
    # compound_file = "compound.pkl"
    # fname = os.path.join(os.path.dirname(__file__), "../data",  compound_file)
    # df.to_pickle(fname)

    #
    # Create and store study dataframe
    #
    # cmd = "SELECT study_id, subst_id, \
    #     normalised_administration_route, normalised_species, normalised_strain, \
    #     exposure_period_days, report_number \
    #     FROM study;"
    cmd = "SELECT study_design.luid AS study_id, subst_id, STANDARDISED_ROUTE AS normalised_administration_route, \
            STANDARDISED_SPECIES as normalised_species, STANDARDISED_STRAIN AS normalised_strain, EXPOSURE_PERIOD AS exposure_period_days, REPORT_NUMBER \
            FROM STUDY_DESIGN \
            JOIN SUBSTANCE_IDS ON STUDY_DESIGN.STRUCTURE_LUID = SUBSTANCE_IDS.LUID"
    
    df = pd.read_sql(cmd, con=conn)
    df.subst_id = df.subst_id.str.upper()
    df.normalised_administration_route = df.normalised_administration_route.str.capitalize()
    df.normalised_species = df.normalised_species.str.capitalize()
    df.normalised_strain = df.normalised_strain.str.capitalize()
    study_file = "study.pkl"
    fname = os.path.join(os.path.dirname(__file__), "../data",  study_file)
    df.to_pickle(fname)

    #
    # Create and store finding dataframe
    #
    # cmd = "SELECT study_id, \
    #             relevance,\
    #             observation_normalised, organ_normalised, normalised_sex, \
    #             dose, grade \
	#             FROM findings_all \
    #                 WHERE source = 'HistopathologicalFinding' \
    #                 AND observation_normalised IS NOT NULL \
    #                 AND organ_normalised IS NOT NULL"
    cmd = "SELECT DISTINCT PARENT_LUID AS study_id, RELEVANCE, \
            STANDARDISED_PATHOLOGY AS observation_normalised, STANDARDISED_ORGAN AS organ_normalised, \
            STANDARDISED_SEX AS normalised_sex, DOSE, grade \
            FROM HISTOPATHOLOGICALFI \
            WHERE STANDARDISED_PATHOLOGY IS NOT NULL \
            AND STANDARDISED_ORGAN IS NOT NULL"
    hp_df = pd.read_sql(cmd, con=conn)
    hp_df['source'] = 'Histopathological'

    cmd = "SELECT DISTINCT PARENT_LUID AS study_id, RELEVANCE, \
            clinical_sign AS observation_normalised, \
            STANDARDISED_SEX AS normalised_sex, DOSE, grade \
            FROM clinical_signs \
            WHERE clinical_sign IS NOT NULL"
    cs_df = pd.read_sql(cmd, con=conn)
    cs_df['source'] = 'ClinicalSigns'

    cmd = "SELECT DISTINCT PARENT_LUID AS study_id, RELEVANCE, \
            stand_gross_pathology AS observation_normalised, \
            STANDARDISED_ORGAN AS organ_normalised, \
            STANDARDISED_SEX AS normalised_sex, DOSE, grade \
            FROM gross_necropsy \
            WHERE stand_gross_pathology IS NOT NULL \
            AND STANDARDISED_ORGAN IS NOT NULL"
    gn_df = pd.read_sql(cmd, con=conn)
    gn_df['source'] = 'GrossNecropsy'

    df = pd.concat([hp_df, cs_df, gn_df])

    df.relevance = df.relevance.str.capitalize()
    df.relevance = df.relevance.fillna('NA')
    df.observation_normalised = df.observation_normalised.str.capitalize()
    df.observation_normalised = df.observation_normalised.fillna('NA')
    df.organ_normalised = df.organ_normalised.str.capitalize()
    df.organ_normalised = df.organ_normalised.fillna('NA')
    df.normalised_sex = df.normalised_sex.str.upper()
    df.grade = df.grade.str.capitalize()
    df.grade = df.grade.fillna('NA')
    find_file = 'findings.pkl.gz'
    fname = os.path.join(os.path.dirname(__file__), "../data",  find_file)
    df.to_pickle(fname, compression='gzip')

    #
    # Create and store anatomy and histopathology ontology dataframe
    #
    cmd = '( WITH RECURSIVE ontology  AS ( \
        select relations."ONTOLOGY_TERM_ID" as child_term, \
        relations."ONTOLOGY_TERM_ID" as parent_term \
        from  input_onto_etox_ontology_relationships relations \
        UNION \
        select onto.child_term as child_term,\
        relations."RELATED_ONTOLOGY_TERM_ID" as parent_term \
        from   input_onto_etox_ontology_relationships relations \
        inner join ontology onto on (onto.parent_term=relations."ONTOLOGY_TERM_ID")\
        )\
        SELECT terms."TERM_NAME" as child_term, \
            terms2."TERM_NAME" as parent_term, \
            terms."ONTOLOGY_NAME" as ontology \
        FROM ontology onto \
        INNER JOIN  input_onto_etox_ontology_terms terms on ( onto.child_term = terms."ONTOLOGY_TERM_ID") \
        INNER JOIN  input_onto_etox_ontology_terms terms2 on ( onto.parent_term= terms2."ONTOLOGY_TERM_ID") \
        WHERE terms."ONTOLOGY_NAME" =\'anatomy\' and  terms2."ONTOLOGY_NAME" =\'anatomy\' \
    )UNION ( \
    WITH RECURSIVE ontology  AS ( \
        select relations."ONTOLOGY_TERM_ID" as child_term, \
            relations."ONTOLOGY_TERM_ID" as parent_term \
            from  input_onto_etox_ontology_relationships relations \
        UNION \
        select onto.child_term as child_term, \
            relations."RELATED_ONTOLOGY_TERM_ID" as parent_term \
            from   input_onto_etox_ontology_relationships relations \
            inner join ontology onto on (onto.parent_term=relations."ONTOLOGY_TERM_ID") \
        ) \
        SELECT terms."TERM_NAME" as child_term, \
            terms2."TERM_NAME" as parent_term, \
            terms."ONTOLOGY_NAME" as ontology \
        FROM ontology onto \
        INNER JOIN  input_onto_etox_ontology_terms terms on ( onto.child_term = terms."ONTOLOGY_TERM_ID") \
        INNER JOIN  input_onto_etox_ontology_terms terms2 on ( onto.parent_term= terms2."ONTOLOGY_TERM_ID") \
        WHERE terms."ONTOLOGY_NAME" =\'histopathology\' and  terms2."ONTOLOGY_NAME" = \'histopathology\')'
           
    df = pd.read_sql(cmd, con=conn)
    df.child_term = df.child_term.str.capitalize()
    df.parent_term = df.parent_term.str.capitalize()
    find_file = "ontology.pkl"
    fname = os.path.join(os.path.dirname(__file__), "../data",  find_file)
    df.to_pickle(fname)

    #
    # Create and store normalization lookup table
    #    
    cmd = 'SELECT syn."VX_VALUE" AS verbatim, \
        term."TERM_NAME" AS normalised \
        FROM public.input_onto_vx_synonyms AS syn \
        JOIN public.input_onto_etox_ontology_terms AS term \
        ON syn."ONTOLOGY_TERM_ID" = term."ONTOLOGY_TERM_ID" \
        WHERE syn."SYNONYM_STATUS" = \'APPROVED\' \
        AND syn."SYNONYM_DOMAIN" IN (\'strain\', \'sex\', \
            \'administration route\', \'species\', \'organ tissue\', \
            \'microscopic finding\');'
    curs.execute(cmd)
    norm_d = {verbatim: normalised for (verbatim, normalised) in curs}

    all_verbatim = list(norm_d.keys())[:]
    for verbatim in all_verbatim:
        verbatim_upper = verbatim.upper()
        norm_d[verbatim_upper] = norm_d[verbatim]
    all_norm = list(set(norm_d.values()))[:]
    for normalised in all_norm:
        normalised_upper = normalised.upper()
        norm_d[normalised_upper] = normalised
    
    norm_file = 'normalisation.pkl'
    fname = os.path.join(os.path.dirname(__file__), '../data',  norm_file)
    with open(fname, 'wb') as f:
        pickle.dump(norm_d, f)


def main ():
    """
    Parse arguments and load the extraction filters.
    """
    parser = argparse.ArgumentParser(description='Create and store \
                compounds, study, findings and ontologies dataframes.')

    parser.add_argument('-u', '--user', help='PostgreSQL DB username.', 
            default='postgres', required=False)
    parser.add_argument('-p', '--passw', help='PostgreSQL DB password.',
            required=True)

    args = parser.parse_args()
    run(args)

if __name__ == '__main__':    
    main()
