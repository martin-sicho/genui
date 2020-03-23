import React from 'react';
import { ComponentWithPagedResources } from '../../../index';

function MoleculeActivityProvider(props) {
  const mol = props.mol;
  const activitySets = props.activitySets;
  const ListComp = props.component;

  const definition = {};
  Object.keys(activitySets).forEach(actSetID => {
    definition[actSetID] = new URL(`${mol.id}/activities/?activity_set=${actSetID}`, props.apiUrls.compoundsRoot);
  });

  return (
    <React.Fragment>
      {/*<h4>Activity Data</h4>*/}
      <ComponentWithPagedResources
        {...props}
        definition={definition}
        // mol={mol}
        // updateCondition={(prevProps, currentProps) => {
        //   return prevProps.mol && (prevProps.mol.id !== currentProps.mol.id)
        // }}
      >
        {
          (activities) => (
            <ListComp {...props} activities={activities}/>
          )
        }
      </ComponentWithPagedResources>
    </React.Fragment>
  )
}

export default MoleculeActivityProvider;