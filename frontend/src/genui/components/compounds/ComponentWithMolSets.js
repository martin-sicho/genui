import React from 'react';

class ComponentWithMolSets extends React.Component {

  constructor(props) {
    super(props);

    this.urlRoots = {
      genericList : new URL('all/', this.props.apiUrls.compoundSetsRoot)
    };

    this.state = {
      fetchUpdates : true,
      compoundSets : null,
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.currentProject && this.state.fetchUpdates) {
      this.fetchUpdates(this.props.currentProject);
    }
  }

  fetchUpdates = (project) => {
    const params = new URLSearchParams();
    params.append('project_id', project.id);
    fetch(this.urlRoots.genericList.toString() + "?" + params.toString())
      .then(response => response.json())
      .then(this.getMolSets)
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
      fetchUpdates : false
    })
  };

  handleAddMolSet = (className, data) => {
    this.setState(prevState => {
      prevState.compoundSets[className].push(data);
      return {
        compoundSets : prevState.compoundSets
      };
    });
  };

  handleAddMolSetList = (className, molsetList, overwrite=false) => {
    this.setState((prevState) => {
      const old_sets = prevState.compoundSets;
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

  handleMolSetDelete = (className, molset) => {
    fetch(this.urlRoots.genericList.toString() + molset.id + '/', {method: 'DELETE'})
      .then(
        () => {
          this.setState(prevState => {
            const molset = prevState.compoundSets[className];
            const idx_del = molset.findIndex(item => item.id === molset.id);
            molset.splice(idx_del, 1);
            prevState.compoundSets[className] = molset; // TODO: check if this is even needed
            return {
              compoundSets : prevState.compoundSets
            };
          });
        }
      ).catch(
      (error) => console.log(error)
    )
    ;
  };

  render() {
    return (
      <React.Fragment>
        {this.props.render(this.state.compoundSets, this.handleAddMolSetList, this.handleAddMolSet, this.handleMolSetDelete)}
      </React.Fragment>
    )
  }
}

export default ComponentWithMolSets;