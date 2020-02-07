import React, {Component} from 'react';
import ChEMBLGrid from './chembl/ChEMBLGrid';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';
import { ComponentWithObjects } from '../../../genui';

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

function MolSetGrid(props) {

  return (
    <div>View for a molset with an undefined class is not yet available.</div>
  )
}

class CompoundsPage extends Component {
  CLASS_TO_COMPONENT = {
    ChEMBLCompounds : ChEMBLGrid,
    MolSet : MolSetGrid
  };

  constructor(props) {
    super(props);

    this.defaultClass = this.props.defaultClass;
    this.classToComponentNoDefault = Object.keys(this.CLASS_TO_COMPONENT).reduce((object, key) => {
      if (key !== this.defaultClass) {
        object[key] = this.CLASS_TO_COMPONENT[key]
      }
      return object
    }, {});
    this.classToComponent = this.CLASS_TO_COMPONENT;
  }

  componentDidMount() {
    this.props.onHeaderChange(
      <HeaderNav
        {...this.props}
        molSetChoices={Object.keys(this.classToComponentNoDefault)}
        onMolSetChoice={this.props.handleAddMolSetList}
      />
    );
  }

  render() {
    const molsets = this.props.compoundSets;

    if (molsets === null) {
      return <div>Loading...</div>
    }

    const molsetsEmpty = Object.keys(molsets).length === 1 && molsets[this.props.defaultClass].length === 0 && molsets.constructor === Object;
    if (molsetsEmpty) {
      return <div><p>There are currently no compound sets. Start by adding one from the actions menu in the top right.</p></div>;
    }

    const classToComponent = this.classToComponentNoDefault;
    return (
      <div className="compound-set-grids">
        {
          Object.keys(molsets).map(MolSetClass => {
            if (classToComponent.hasOwnProperty(MolSetClass)) {
              const MolsetComponent = classToComponent[MolSetClass];
              return (
                <div key={MolSetClass} className={MolSetClass}>
                  <MolsetComponent
                    {...this.props}
                    molsets={molsets[MolSetClass]}
                    currentMolsetClass={MolSetClass}
                  />
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

class Compounds extends React.Component {
  DEFAULT_CLASS = "MolSet";

  render() {
    const defaultClass = this.DEFAULT_CLASS;
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
              compoundSets={compoundSets}
              defaultClass={defaultClass}
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
