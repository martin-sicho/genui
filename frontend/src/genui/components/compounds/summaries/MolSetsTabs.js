import { TabWidget } from '../../../index';
import React from 'react';
import CompoundList from './CompoundList';

export default function MolSetsTabs(props) {
  const groupedMols = props.groupedMols;
  const molsets = props.molsets;
  const mols = props.mols;

  const tabs = [];
  let activeTab = undefined;
  Object.keys(groupedMols).forEach(molsetID => {
    const molset = molsets.find(item => item.id === Number(molsetID));
    const data = groupedMols[molsetID];
    if (!activeTab) {
      activeTab = molset.name;
    }
    tabs.push({
      title: molset.name
      , renderedComponent : props => (
        <CompoundList
          {...props}
          paginate={true}
          molset={molset}
          mols={data.map(item => item.mol)}
          // points={data.map(item => item.point)}
        />
      )
    })
  });

  return (
    <div className="genui-map-molsets-grouped">
      {
        mols.length > 0 ? <TabWidget {...props} tabs={tabs} activeTab={activeTab}/> : <p>Select molecules in the map to see details.</p>
      }
    </div>
  )
}