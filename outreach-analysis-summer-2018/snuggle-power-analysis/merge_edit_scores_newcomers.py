
# coding: utf-8

# # Prototype: Merging Newcomer Dataframes with ORES scores for Newcomer Contributions
# June 8, 2018 J. Nathan Matias
#
# Using data sources from http://paws-public.wmflabs.org/paws-public/User:Juliakamin/Querying%20new%20editors%20via%20sql.ipynb



import os, time, datetime, csv, glob, math, datetime, pprint
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from dateutil import parser
import click

@click.command()
@click.option('--lang', default='es', help='the wiki language to target')
@click.option('--datadir', default='data', help='where to look for the data')
@click.option('--resultsdir', default='results', help='where to plop results when done')
def main(lang, datadir, resultsdir):
    newcomer_file = os.path.join(datadir, lang+"_newcomer_list.csv")
    newcomer_revisions_files = glob.glob(
        os.path.join(datadir, lang + "_newcomer_revisions*.csv"))



    newcomers = {}
    counter = 0
    with open(newcomer_file, "r") as f:
        for newcomer in csv.DictReader(f.readlines()):
            ## REMOVE OUT THE PANDAS SEQUENTIAL INDEX
            ## IF IT EXISTS
            if '' in newcomer.keys():
                del newcomer['']
            newcomer['wiki']  = lang
            newcomer['registration.date'] = datetime.datetime.strptime(
            newcomer['registration'].replace("b'","").replace("'",""),
            "%Y%m%d%H%M%S")
            newcomers[newcomer['user id']] = newcomer
            counter += 1




    newcomer_revisions = defaultdict(list)
    all_ids = set()
    for filename in newcomer_revisions_files:
        with open(filename, "r") as f:
            print(filename)
            for revision in csv.DictReader(f.readlines()):
                revision['wiki'] = lang
                revision_id = revision['revision id']
                if('' in revision.keys()):
                    del revision['']
                if(revision_id not in all_ids):
                    newcomer_revisions[revision['user id']].append(revision)
                    all_ids.add(revision_id)

    for key, revisions in newcomer_revisions.items():
        newcomer_revisions[key] = sorted(revisions,
                                         key=lambda x:parser.parse(x['revision time']))


    # ### Data Validation
    # Here, we confirm that every newcomer has at least one revision
    # And that there aren't any revisions that have no newcomer



    print("{0} total unique newcomers".format(len(set([x for x in newcomers.keys()]))))
    print("{0} total newcomers with edit records in the dataset".format(
            len(newcomer_revisions.keys())))




    newcomers_in_revision_set = set()
    revisions_not_in_newcomer_set = set()
    newcomers_not_in_revision_set = set()
    for user_id, revisions in newcomer_revisions.items():
        if user_id in newcomers.keys():
            newcomers_in_revision_set.add(user_id)
        else:
            revisions_not_in_newcomer_set.add(user_id)

    for user_id in newcomers.keys():
        if user_id not in newcomer_revisions.keys():
            newcomers_not_in_revision_set.add(key)

    print("{0} newcomers in revision set".format(len(newcomers_in_revision_set)))
    print("{0} revisions not in newcomer set".format(len(revisions_not_in_newcomer_set)))
    print("{0} newcomers not in revision set".format(len(newcomers_not_in_revision_set)))


    # ### Create Revision Dataframe with Information on users


    all_revisions = []
    for user_id, revisions in newcomer_revisions.items():
        newcomer = newcomers[user_id]
        first_revision_time = None
        if(len(revisions)>0):
            first_revision_time = parser.parse(revisions[0]['revision time'])
            treat_time = first_revision_time + datetime.timedelta(hours=48)
        for revision in revisions:
            revision_time = parser.parse(revision['revision time'])
            revision['registration'] = newcomer['registration.date']
            revision['days.since.registration'] = (revision_time -
                                                   revision['registration']).days
            revision['days.since.simulated.treat'] = (revision_time - treat_time).days
            revision['edits.6.months'] = newcomer['edit count']
            all_revisions.append(revision)

    all_revs_file = os.path.join(resultsdir, lang+"_revisions_with_user_11.2017.csv")
    pd.DataFrame(all_revisions).to_csv(all_revs_file)
    print('saved file: {}'.format(all_revs_file))


    # ### Create User Dataframe with Summary Stats on Revisions
    # For the power analysis, assume that the participant will be 'treated' within 2 days (48 hrs) of making their first edit. Then count the following information starting 2 days after their first edit:
    # * edits in the 7 day period between 3-4 weeks
    # * edits in the 4 week period between 1-2 months
    # * edits in the 12 week period betwen 2-6 months


    for user_id, newcomer in newcomers.items():
        for key in ['first.n.goodfaith.mean',
                    'first.n.goodfaith.median',
                    'first.n.damaging.mean',
                    'first.n.damaging.median',
                    'first.goodfaith.score',
                    'first.damaging.score']:
            newcomer[key] = None

        newcomer['ORES.rated.revisions'] = 0

        ## edits counted from registration time
        newcomer['edits.6.months'] = 0
        newcomer['edits.2.weeks']  = 0
        newcomer['edits.4.weeks'] = 0
        newcomer['edits.8.weeks'] = 0
        newcomer['edits.12.weeks'] = 0

        ## time interval simulation for power analysis
        ## following the outcome variables used in this paper:
        ## https://osf.io/preprints/socarxiv/8qsv6/
        newcomer['edits.3.4.weeks'] = 0
        newcomer['edits.4.8.weeks'] = 0
        newcomer['edits.8.24.weeks'] = 0

        ## SURVIVAL MEASURE OVER 20 WEEK PERIODS (5 MONTHS)
        ## BASED ON WHETHER THEY MADE AT LEAST ONE EDIT
        ## AT ANY TIME AFTER THE OBSERVED PERIOD
        ## This measure uses days.since.simulated.treat
        ## for calculation
        ## We only consider 5 months, since we only have six
        ## months of information, and since some joined at the end of Nov
        survival_weeks = 20
        for i in range(1, survival_weeks+1):
            newcomer['survival.week.period.' + str(i)] = False

        if user_id in newcomer_revisions:
            revisions = newcomer_revisions[user_id]

            ### COUNT REVISIONS OVER A PERIOD OF TIME
            newcomer['edits.6.months'] = len(revisions)
            newcomer['edits.2.weeks'] = len([x for x in revisions if x['days.since.registration'] <= 7*2])
            newcomer['edits.4.weeks'] = len([x for x in revisions if x['days.since.registration'] <= 7*4])
            newcomer['edits.8.weeks'] = len([x for x in revisions if x['days.since.registration'] <= 7*8])
            newcomer['edits.12.weeks'] = len([x for x in revisions if x['days.since.registration'] <= 7*12])

            ## SET TIME INTERVAL SIMULATION FOR POWER ANALYSIS
            newcomer['edits.3.4.weeks'] = len([x for x in revisions if
                                               x['days.since.simulated.treat'] >= 7*3 and
                                               x['days.since.simulated.treat'] < 7*4 ])
            newcomer['edits.4.8.weeks'] = len([x for x in revisions if
                                               x['days.since.simulated.treat'] >= 7*4 and
                                               x['days.since.simulated.treat'] < 7*8 ])
            newcomer['edits.8.24.weeks'] = len([x for x in revisions if
                                               x['days.since.simulated.treat'] >= 8*4 and
                                               x['days.since.simulated.treat'] < 7*24 ])

            ## SET SURVIVAL COLUMNS
            ## first, find eligible revisions
            ## Negative day scores are for edits made in the 48 hours from the point of the first edit
            ## and include the first edit
            survival_revision_days = [x['days.since.simulated.treat'] for x in revisions if
                                  math.ceil(x['days.since.simulated.treat'] / 7) <= survival_weeks and
                                  x['days.since.simulated.treat'] > 0]
            ## If there is at least one eligible revision
            ## then we should update all appropriate survival periods to True
            if(len(survival_revision_days)>0):
                final_revision_day = survival_revision_days[-1]
                final_revision_period = math.ceil(final_revision_day / 7)
                for i in range(1, final_revision_period+1):
                    newcomer['survival.week.period.' + str(i)] = True

            ### AGGREGATE ORES SCORES
            eligible_revisions = [x for x in revisions if
                                  x['goodfaith']!='na' and x['damaging']!='na' and
                                  x['goodfaith']!='error' and x['damaging']!='error']
            newcomer['ORES.rated.revisions'] = len(eligible_revisions)
            if(len(eligible_revisions)>0):
                newcomer['first.goodfaith.score'] = float(eligible_revisions[0]['goodfaith'])
                newcomer['first.damaging.score'] = float(eligible_revisions[0]['damaging'])
                newcomer['first.n.goodfaith.mean'] = np.mean([float(x['goodfaith']) for x in eligible_revisions])
                newcomer['first.n.goodfaith.median'] = np.median([float(x['goodfaith']) for x in eligible_revisions])
                newcomer['first.n.damaging.mean'] = np.mean([float(x['damaging']) for x in eligible_revisions])
                newcomer['first.n.damaging.median'] = np.median([float(x['damaging']) for x in eligible_revisions])



    ores_file = os.path.join(resultsdir, lang+"_newcomers_with_ores_scores_11.2017.csv")
    pd.DataFrame(list(newcomers.values())).to_csv(ores_file)
    print('saved file: {}'.format(ores_file))


    # ### Make Survival Dataframe

    survival_week_records = []
    for user_id, newcomer in newcomers.items():
        if(int(newcomer['edit count'])>0):
            for i in range(1, survival_weeks+1):
                survival_week_records.append({"user.id": user_id,
                                       "week":i,
                                       "first.goodfaith.score":newcomer['first.goodfaith.score'],
                                       "first.damaging.score":newcomer['first.damaging.score'],
                                       "survived":int(newcomer["survival.week.period."+str(i)])})
    print("Generated {0} survival week periods".format(len(survival_week_records)))

    survival_file = os.path.join(resultsdir, lang+"_newcomer_survival_week_periods.2017.csv")
    pd.DataFrame(survival_week_records).to_csv(survival_file)
    print('saved file: {}'.format(survival_file))

if __name__ == '__main__':
    main()
    print('Finshed. ☺')
