import React from 'react';

export function ActivitiesList(props) {
  const set = props.set;
  const activities = props.activities;

  return (
    <React.Fragment>
      <h6>{set.name}</h6>
      <ul>
        {
          activities.map(item => {
            return (
              <li key={item.id}>{item.value}</li>
            )
          })
        }
      </ul>
    </React.Fragment>
  )
}

export function ActivitySetList(props) {
  const actSets = props.activitySets;
  const activities = props.activities;

  return (
    <ul>
      {
        Object.keys(actSets).map(key => {
          const set = actSets[key];
          if (activities[key].length > 0) {
            return <li key={set.id}><ActivitiesList set={set} activities={activities[key]}/></li>
          } else {
            return null;
          }
        })
      }
    </ul>
  )
}