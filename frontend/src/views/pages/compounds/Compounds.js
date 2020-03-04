import React from 'react';
import { ComponentWithObjects, CompoundsPage } from '../../../genui';
import ChEMBLGrid from './chembl/ChEMBLGrid';
import GeneratedGrid from './generated/GeneratedGrid';

class Compounds extends React.Component {
  CLASS_TO_COMPONENT = {
    ChEMBLCompounds : ChEMBLGrid,
    GeneratedMolSet : GeneratedGrid,
  };

  render() {
    const defaultClass = "MolSet";
    return (
      <ComponentWithObjects
        {...this.props}
        objectListURL={new URL('all/', this.props.apiUrls.compoundSetsRoot)}
        emptyClassName={defaultClass}
        render={
          (
            compoundSets,
            handleAddMolSetList,
            handleAddMolSet,
            handleMolSetDelete,
          ) => {
            return (<CompoundsPage
              {...this.props}
              classToComponentMap={this.CLASS_TO_COMPONENT}
              compoundSets={compoundSets}
              defaultClass={defaultClass}
              ignoreDefault={true}
              handleAddMolSetList={handleAddMolSetList}
              handleAddMolSet={handleAddMolSet}
              handleMolSetDelete={handleMolSetDelete}
            />)
          }
        }
      />
    )
  }
}

export default Compounds;
