import React from 'react';
import { Table } from 'reactstrap';
import { ComponentWithResources } from '../../../index';

export function ActivitySetStatsTable(props) {

  const [clickedID, setClickedID] = React.useState(null);
  return (
    <Table size="sm" responsive hover>
      <thead>
      <tr>
        <th>Activity Type</th>
        <th>Data Points</th>
        <th>Molecules</th>
        <th>Activity Set</th>
      </tr>
      </thead>
      <tbody>
      {
        props.summaries.map(summary => {
          const Data = () => (
            <React.Fragment>
              <td>{summary.type.value}</td>
              <td>{summary.activities}</td>
              <td>{summary.molecules}</td>
              <td>{summary.activitySet.name}</td>
            </React.Fragment>
          );
          if (props.selectable) {
            return (
              <tr key={summary.id} className={clickedID && (clickedID === summary.id) ? "bg-success text-dark" : null} onClick={() => {setClickedID(summary.id); props.onSelect(summary)}}>
                <Data/>
              </tr>
            )
          } else {
            return (
              <tr key={summary.id}>
                <Data/>
              </tr>
            )
          }
        })
      }
      </tbody>
    </Table>
  )
}

export function MolsetActivitiesSummary(props) {
  const definition_summary = {};
  props.molset.activities.forEach(actsetid => definition_summary[actsetid] = new URL(`${actsetid}/summary/`, props.apiUrls.activitySetsRoot));

  const updateCond = (prevProps, nextProps) => prevProps.molset.id !== nextProps.molset.id;
  return (
    <ComponentWithResources
      {...props}
      definition={definition_summary}
      updateCondition={updateCond}
    >
      {
        (summaryLoaded, summaries) => {
          const definition_actsets = {};
          props.molset.activities.forEach(actsetid => definition_actsets[actsetid] = new URL(`${actsetid}/`, props.apiUrls.activitySetsRoot));
          return (
            <ComponentWithResources
              {...props}
              definition={definition_actsets}
              updateCondition={updateCond}
            >
              {(actsetsLoaded, actsets) => {
                if (summaryLoaded && actsetsLoaded) {
                  const items = [];
                  Object.keys(summaries).forEach(actsetID => {
                    const actset = actsets[actsetID];
                    summaries[actsetID].typeSummaries.forEach(summary => {
                      items.push({
                        id: `${actset.id}_${summary.type.value}`,
                        name: `${summary.type.value} from ${actset.name}`,
                        activitySet: actset,
                        type: summary.type,
                        molecules: summary.moleculesTotal,
                        activities: summary.activitiesTotal
                      })
                    })
                  });
                  return (
                    <React.Fragment>
                      {props.message ? <p>{props.message}</p> : null}
                      <ActivitySetStatsTable
                        {...props}
                        summaries={items}
                        activitySets={actsets}
                      />
                    </React.Fragment>
                  )
                } else {
                  return <p>Fetching data...</p>
                }
              }}
            </ComponentWithResources>
          )
        }
      }
    </ComponentWithResources>
  )
}