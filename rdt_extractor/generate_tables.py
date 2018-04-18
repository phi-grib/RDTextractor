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

def run(args):

    conn = psycopg2.connect(host='172.20.16.76', dbname='vitic2016_1', 
                            user=args.user, password=args.passw)
    curs = conn.cursor()

    #
    # Create and store study dataframe
    #
    cmd = "SELECT study_id, subst_id, normalised_sex, \
        normalised_administration_route, normalised_species, \
        exposure_period_days, report_number \
        FROM study;"
    
    df = pd.read_sql(cmd, con=conn)
    study_file = "study.pkl"
    fname = os.path.join(os.path.dirname(__file__), "../data", study_file)
    df.to_pickle(fname)

    #
    # Create and store finding dataframe
    #
    cmd = "SELECT study_id, \
                (CASE \
                WHEN relevance IS NULL THEN 'NA' \
                WHEN relevance = 'Treatment related' THEN 'treatment related' \
                ELSE relevance END) AS relevance,\
                observation_normalised, organ_normalised, dose \
	            FROM findings_all \
                    WHERE source = 'HistopathologicalFinding' \
                    AND observation_normalised IS NOT NULL \
                    AND organ_normalised IS NOT NULL"
    df = pd.read_sql(cmd, con=conn)
    df.relevance = df.relevance.fillna('NA')
    df.observation_normalised = df.observation_normalised.fillna('NA')
    df.organ_normalised = df.organ_normalised.fillna('NA')
    find_file = 'findings.pkl.gz'
    fname = os.path.join(os.path.dirname(__file__), "../data", find_file)
    df.to_pickle(fname, compression='gzip')

    #
    # Create and store expanded anatomy ontology dataframe
    #
    cmd = '( WITH RECURSIVE ontology AS ( \
        select relations."ONTOLOGY_TERM_ID" as child_term, \
        relations."ONTOLOGY_TERM_ID" as parent_term \
        from input_onto_etox_ontology_relationships relations \
        UNION \
        select onto.child_term as child_term,\
        relations."RELATED_ONTOLOGY_TERM_ID" as parent_term \
        from input_onto_etox_ontology_relationships relations \
        inner join ontology onto on (onto.parent_term=relations."ONTOLOGY_TERM_ID")\
        )\
        SELECT terms."TERM_NAME" as child_term, \
            terms2."TERM_NAME" as parent_term, \
            terms."ONTOLOGY_NAME" as ontology \
        FROM ontology onto \
        INNER JOIN input_onto_etox_ontology_terms terms on ( onto.child_term = terms."ONTOLOGY_TERM_ID") \
        INNER JOIN input_onto_etox_ontology_terms terms2 on ( onto.parent_term= terms2."ONTOLOGY_TERM_ID") \
        WHERE terms."ONTOLOGY_NAME" =\'anatomy\' and terms2."ONTOLOGY_NAME" =\'anatomy\')'
           
    df = pd.read_sql(cmd, con=conn)
    find_file = "anatomy_ontology.pkl"
    fname = os.path.join(os.path.dirname(__file__), "../data", find_file)
    df.to_pickle(fname)

    #
    # Create and store anatomy and histopathology ontology dataframe
    #
    cmd = '(WITH RECURSIVE ontology AS ( \
        select relations."ONTOLOGY_TERM_ID" as child_term, \
            relations."ONTOLOGY_TERM_ID" as parent_term \
            from input_onto_etox_ontology_relationships relations \
        UNION \
        select onto.child_term as child_term, \
            relations."RELATED_ONTOLOGY_TERM_ID" as parent_term \
            from input_onto_etox_ontology_relationships relations \
            inner join ontology onto on (onto.parent_term=relations."ONTOLOGY_TERM_ID") \
        ) \
        SELECT terms."TERM_NAME" as child_term, \
            terms2."TERM_NAME" as parent_term, \
            terms."ONTOLOGY_NAME" as ontology \
        FROM ontology onto \
        INNER JOIN input_onto_etox_ontology_terms terms on ( onto.child_term = terms."ONTOLOGY_TERM_ID") \
        INNER JOIN input_onto_etox_ontology_terms terms2 on ( onto.parent_term= terms2."ONTOLOGY_TERM_ID") \
        WHERE terms."ONTOLOGY_NAME" =\'histopathology\' and terms2."ONTOLOGY_NAME" = \'histopathology\')'
           
    df = pd.read_sql(cmd, con=conn)
    df = df[df.parent_term != 'morphologic change']
    find_file = "morph_changes_ontology.pkl"
    fname = os.path.join(os.path.dirname(__file__), "../data", find_file)
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
    fname = os.path.join(os.path.dirname(__file__), '../data', norm_file)
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
