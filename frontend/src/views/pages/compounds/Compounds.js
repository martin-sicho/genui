import React, {Component} from 'react';

class Compounds extends Component {

  constructor(props) {
    super(props);

    this.urlRoots = {
      genericList : new URL('generic/', this.props.apiUrls.compoundSetsRoot),
      ChEMBLCompounds : new URL('chembl/', this.props.apiUrls.compoundSetsRoot),
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
    this.setState({
      compoundSets : data,
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

    if (!molsets) {
      return <div>Loading...</div>
    }

    return (
      <div>
        <span>First compound set: {molsets[0].name}</span>
      </div>
    );
  }
}

export default Compounds;
