import React from 'react';
import './compound-list-styles.css'
import { resolve } from '../../../utils';

export function DataPair(props) {
  const DataComponent = props.dataComponent;
  return DataComponent ? (
    <p><b>{props.title}:</b> <DataComponent {...props}/></p>
  ) : (
    <p><b>{props.title}:</b> {props.data}</p>
  )
}

export function MoleculeMetadata(props) {
  const mol = props.mol;

  const extraData = [];
  if (props.extraInfoFields) {
    props.extraInfoFields.forEach(definition => {
      if (mol.className === definition.className) {
        // TODO: check if data items and prop names have the same length
        const data = definition.dataItems.map(dataItem => resolve(dataItem, mol));
        const sentProps = {};
        definition.propNames.forEach((propName, index) => sentProps[propName] = data[index]);
        extraData.push({
          title: definition.displayName,
          component: definition.component,
          props: sentProps,
          className: definition.className,
        })
      }
    });
  }

  return (
    <React.Fragment>
      <DataPair title="SMILES" data={mol.smiles}/>
      <DataPair title="InChiKey" data={mol.inchiKey}/>
      {
        extraData.map(data => (
          <DataPair key={data.title} title={data.title} dataComponent={data.component} {...data.props} />
        ))
      }
    </React.Fragment>
  )
}