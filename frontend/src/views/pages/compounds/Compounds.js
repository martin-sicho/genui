import React, {Component} from 'react';
import ChEMBLGrid from './chembl/ChEMBLGrid';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';
import { ComponentWithMolSets } from '../../../genui';

function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Actions
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Add New...</DropdownToggle>
        <DropdownMenu>
          {
            props.molSetChoices.map(choice =>
              (<DropdownItem
                key={choice}
                onClick={() => {props.onMolSetChoice(choice, [])}}
              >
                {choice}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
    </DropdownMenu>
  </UncontrolledDropdown>)
}

class CompoundsPage extends Component {

  CLASS_TO_COMPONENT = {
    ChEMBLCompounds : ChEMBLGrid
  };

  componentDidMount() {
    this.props.onHeaderChange(<HeaderNav {...this.props} molSetChoices={Object.keys(this.CLASS_TO_COMPONENT)} onMolSetChoice={this.props.handleAddMolSetList}/>);
  }

  render() {
    const molsets = this.props.compoundSets;

    if (molsets === null) {
      return <div>Loading...</div>
    }

    if (Object.keys(molsets).length === 0) {
      return <p>Start by selecting a compound set to create from the action menu.</p>
    }

    return (
      <div className="compound-set-grids">
        {
          Object.keys(molsets).map(MolSetClass => {
            if (this.CLASS_TO_COMPONENT.hasOwnProperty(MolSetClass)) {
              const MolsetComponent = this.CLASS_TO_COMPONENT[MolSetClass];
              return (
                <div key={MolSetClass} className={MolSetClass}>
                  <MolsetComponent {...this.props} molsets={molsets[MolSetClass]} currentMolsetClass={MolSetClass}/>
                </div>
              )
            } else {
              console.log(`Unknown class: ${MolSetClass}`);
              return null;
            }
          })
        }
      </div>
    );
  }
}

function Compounds(props) {
  return (
    <ComponentWithMolSets
      {...props}
      render={
        (
          compoundSets,
          handleAddMolSetList,
          handleAddMolSet,
          handleMolSetDelete,
        ) => {
          return (<CompoundsPage
            {...props}
            compoundSets={compoundSets}
            handleAddMolSetList={handleAddMolSetList}
            handleAddMolSet={handleAddMolSet}
            handleMolSetDelete={handleMolSetDelete}
          />)
        }
      }
    />
  )
}

export default Compounds;
