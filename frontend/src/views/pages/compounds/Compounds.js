import React, {Component} from 'react';
import ChEMBLGrid from './chembl/ChEMBLGrid';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

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

class Compounds extends Component {

  CLASS_TO_COMPONENT = {
    ChEMBLCompounds : ChEMBLGrid
  };

  // FIXME: this page is empty when there are no existing compound sets to display
  // add buttons to the toolbar to initialize new compund set groups with the proper forms

  constructor(props) {
    super(props);

    this.urlRoots = {
      genericList : new URL('all/', this.props.apiUrls.compoundSetsRoot)
    };

    this.state = {
      isLoading : true,
      fetchUpdates : true,
      compoundSets : null,
    }
  }

  // custom methods

  fetchUpdates = () => {
    if (this.props.currentProject && this.state.fetchUpdates) {
      const project = this.props.currentProject;
      const params = new URLSearchParams();
      params.append('project_id', project.id);
      fetch(this.urlRoots.genericList.toString() + "?" + params.toString())
        .then(response => response.json())
        .then(this.getMolSets)
    }
  };

  getMolSets = (data) => {
    const compoundSets = {};
    for (const cset of data) {
      if (!compoundSets.hasOwnProperty(cset.className)) {
        compoundSets[cset.className] = [];
      }
      compoundSets[cset.className].push(cset);
    }
    this.setState({
      compoundSets : compoundSets,
      fetchUpdates : false,
      isLoading : false
    })
  };

  handleMolSetListAdd = (className, molsetList, overwrite=false) => {
    this.setState((prev) => {
      const old_sets = prev.compoundSets;
      if (old_sets.hasOwnProperty(className)) {
        if (overwrite) {
          old_sets[className] = molsetList;
        } else {
          old_sets[className] = old_sets[className].concat(molsetList);
        }
      } else {
        old_sets[className] = molsetList;
      }
      return {
        compoundSets : old_sets
      }
    });
  };

  // component methods

  componentDidMount() {
    this.props.onHeaderChange(<HeaderNav {...this.props} molSetChoices={Object.keys(this.CLASS_TO_COMPONENT)} onMolSetChoice={this.handleMolSetListAdd}/>);
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
      this.fetchUpdates();
  }

  render() {
    const molsets = this.state.compoundSets;

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
                  <MolsetComponent {...this.props} molsets={molsets[MolSetClass]}/>
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

export default Compounds;
