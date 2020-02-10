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

class CompoundsPage extends Component {

  constructor(props) {
    super(props);

    this.defaultClass = this.props.defaultClass;
    this.ignoreDefault = this.props.ignoreDefault;
    this.classToComponent = this.props.classToComponentMap;
    this.classToComponentNoIgnore = this.ignoreDefault ? Object.keys(this.classToComponent).reduce((object, key) => {
      if (key !== this.defaultClass) {
        object[key] = this.classToComponent[key]
      }
      return object
    }, {}) : undefined;
  }

  componentDidMount() {
    this.props.onHeaderChange(
      <HeaderNav
        {...this.props}
        molSetChoices={Object.keys(this.ignoreDefault ? this.classToComponentNoIgnore : this.classToComponent)}
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

    const classToComponent = this.classToComponentNoIgnore ? this.classToComponentNoIgnore : this.classToComponent;
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
              console.log(`Ignored class without a component: ${MolSetClass}`);
              return null;
            }
          })
        }
      </div>
    );
  }
}

class Compounds extends React.Component {
  CLASS_TO_COMPONENT = {
    ChEMBLCompounds : ChEMBLGrid
    // GeneratedMolSet : GeneratedMolSetGrid,
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
