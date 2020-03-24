import React from 'react';
import { ComponentWithResources } from '../../../index';

export default function MoleculePropsProvider(props) {
  const mol = props.mol;
  const Component = props.component;
  const propsList = props.propsList;

  return (
    <React.Fragment>
      {/*<h4>Activity Data</h4>*/}
      <ComponentWithResources
        {...props}
        definition={{
          mol: new URL(`${mol.id}/?properties=${propsList.join(',')}`, props.apiUrls.compoundsRoot)
        }}
      >
        {
          (complete, molWithProps) => (
            complete ? <Component {...props} molWithProperties={molWithProps.mol}/> : null
          )
        }
      </ComponentWithResources>
    </React.Fragment>
  )
}