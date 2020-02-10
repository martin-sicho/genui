import React from 'react';
import ChEMBLInfo from './tabs/ChEMBLInfo';
import ChEMBLCompounds from './tabs/ChEMBLCompounds';
import {GenericMolSetCard} from '../../../../genui';

function ChEMBLCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : ChEMBLInfo,
    },
    {
      title: "Molecules",
      renderedComponent : ChEMBLCompounds
    }
  ];

  const url = new URL(`chembl/${props.molset.id}/`, props.apiUrls.compoundSetsRoot);
  return (
    <GenericMolSetCard {...props} tabs={tabs} molsetURL={url}/>
  )
}

export default ChEMBLCard;