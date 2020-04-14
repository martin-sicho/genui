import { ComponentWithPagedResources, MolsToMolSetGroups } from '../../../../genui';
import React from 'react';

export default function ActivitiesAggregator(props) {
  return (
    <MolsToMolSetGroups
      {...props}
    >
      {
        (groups) => {
          const groupedMols = groups;
          const definition = {};
          // const molsetToActivities = {};
          Object.keys(groupedMols).forEach(molsetID => {
            const molset = props.molsets.find(molset => molset.id === Number(molsetID));
            const mols =  groupedMols[molsetID].map(mol => mol.mol);
            // molsetToActivities[molsetID] = {};
            // molsetToActivities[molsetID]['molecules'] = mols;
            // molsetToActivities[molsetID]['molset'] = molset;
            // molsetToActivities[molsetID]['acsets'] = {};
            Object.keys(props.activitySets).forEach(key => {
              if (molset.activities.includes(Number(key))) {
                // molsetToActivities[molsetID]['acsets'][key] = {};
                // molsetToActivities[molsetID]['acsets'][key]['acset'] = props.activitySets[key];
                // molsetToActivities[molsetID]['acsets'][key]['activities'] = [];
                definition[key] = new URL(`${key}/activities/?mols=${mols.map(mol => mol.id).join(',')}`, props.apiUrls.activitySetsRoot);
              }
            });
          });

          // console.log(definition);
          return (
            <ComponentWithPagedResources
              {...props}
              definition={definition}
            >
              {
                (data, finished, revision) => {
                  if (finished) {
                    // console.log(finished);
                    // console.log(revision);
                    // console.log(data);
                    // console.log('----');

                    const activitiesAll = [];
                    Object.keys(data).forEach(key => {
                      data[key].forEach(item => activitiesAll.push(item))
                    });

                    // Object.keys(groupedMols).forEach(molsetID => {
                    //   Object.keys(molsetToActivities[molsetID]['acsets']).forEach(acsetID => {
                    //     let activities = data[acsetID];
                    //     activities = groupBy(activities, 'type.id');
                    //     molsetToActivities[molsetID]['acsets'][acsetID]['activities'] = activities;
                    //   });
                    // });

                    // console.log(activitiesAll);

                    return props.children(activitiesAll);
                  } else {
                    return props.children(null)
                  }
                }
              }
            </ComponentWithPagedResources>
          )
        }
      }
    </MolsToMolSetGroups>
  )
}