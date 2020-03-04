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

  return (
    <GenericMolSetCard {...props} tabs={tabs}/>
  )
}

export default ChEMBLCard;