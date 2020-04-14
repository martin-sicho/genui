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

    let extraInfoFields = [];
    if (props.molsetClassToExtraInfoFields && props.molsetClassToExtraInfoFields.hasOwnProperty(molset.className)) {
      extraInfoFields = props.molsetClassToExtraInfoFields[molset.className];
    }

    let extraActivityFields = [];
    if (props.molsetClassToExtraActivityFields && props.molsetClassToExtraActivityFields.hasOwnProperty(molset.className)) {
      extraActivityFields = props.molsetClassToExtraActivityFields[molset.className];
    }

    let molsetListUrl = null;
    if (props.molsetClassToURLs && props.molsetClassToURLs.hasOwnProperty(molset.className)) {
      molsetListUrl = props.molsetClassToURLs[molset.className];
    }

    tabs.push({
      title: molset.name
      , renderedComponent : props => (
        <CompoundList
          {...props}
          paginate={true}
          molset={molset}
          mols={data.map(item => item.mol)}
          extraInfoFields={extraInfoFields}
          extraActivityFields={extraActivityFields}
          molsetListUrl={molsetListUrl}
          // points={data.map(item => item.point)}
        />
      )
    })
  });

  return (
    <div className="genui-map-molsets-grouped">
      {
        mols.length > 0 ? <TabWidget {...props} tabs={tabs} activeTab={activeTab}/> : null
      }
    </div>
  )
}