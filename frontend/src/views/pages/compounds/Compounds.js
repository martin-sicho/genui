import React, {Component} from 'react';
import ChEMBLGrid from './chembl/ChEMBLGrid';

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
      compoundSets : null
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

  // component methods

  componentDidUpdate(prevProps, prevState, snapshot) {
      this.fetchUpdates();
  }

  render() {
    const molsets = this.state.compoundSets;

    if (molsets === null) {
      return <div>Loading...</div>
    }

    return (
      <div className="compound-set-grids">
        {Object.keys(molsets).map(MolSetClass => {
          const MolsetComponent = this.CLASS_TO_COMPONENT[MolSetClass];
          return (<div key={MolSetClass} className={MolSetClass}>
            <MolsetComponent {...this.props} molsets={molsets[MolSetClass]}/>
          </div>)
        })}
      </div>
    );
  }
}

export default Compounds;
